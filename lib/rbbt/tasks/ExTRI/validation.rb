module ExTRI


  #{{{ Threshold

  input :score, :float, "Score threshold", 1.6
  input :pmids, :integer, "Min number of PMIDS", 1000
  input :sentences, :integer, "Min number of sentences", 2
  dep :ExTRI_counts, :compute => :produce
  task :threshold => :tsv do |score, pmids, sentences|
    tsv = step(:ExTRI_counts).load
    selected = []
    selected.concat tsv.select("Interaction score"){|s| s.to_f >= score }.keys
    selected.concat tsv.select("PMID counts"){|s| s.to_f >= pmids }.keys
    selected.concat tsv.select("Sentence counts"){|s| s.to_f >= sentences }.keys

    tsv.select(selected)
  end

  dep :threshold
  task :threshold_pairs => :array do 
    step(:threshold).load.keys.collect{|k| k.split(":").values_at(2,3) * ":" }.uniq
  end

  #{{{ Validation datasets

  dep :flagged_tfs
  task :validation_dataset => :tsv do
    flagged_tfs = step(:flagged_tfs).load

    validation_score = Rbbt.data["Astrid Validation"]["validation_samples_260317.tsv"].tsv :type => :list
    validation_lowscore = Rbbt.data["Astrid Validation"]["validation_sample_lowscore_290317.tsv"].tsv :type => :list
    validation_pmid = Rbbt.data["Astrid Validation"]["validation_sample_pmids_280317_AL_124-labels.tsv"].tsv :type => :list

    validation = validation_score

    validation.add_field "NOT TRUE" do
      nil
    end

    validation_lowscore.through do |k,v|
      nv = []
      validation.fields.each_with_index do |f,i|
        begin
          nv[i] = v[f]
        rescue
        end
      end
      validation[k] = nv
    end

    validation_pmid.through do |k,v|
      nv = []
      validation.fields.each_with_index do |f,i|
        begin
          nv[i] = v[f]
        rescue
        end
      end
      validation[k] = nv
    end

    validation.add_field "Set" do |k,v|
      if validation_pmid[k]
        "PMID"
      elsif validation_lowscore[k]
        "Lowscore"
      else
        "Score"
      end
    end

    validation = validation.select{|k,v| tf = k.split(":")[2]; ! flagged_tfs.include? tf }

    validation.add_field "Valid" do |k,v|
      if (v["0"].nil? || v["0"] == "") && (v["NOT TRUE"].nil? || v["NOT TRUE"] == "")
        "Valid"
      else
        "Not Valid"
      end
    end

    validation.key_field = "PMID:Sentence ID:TF:TG"

    validation.fields = validation.fields.collect{|f| f.sub(/(..) Associated Gene Name/, '\1 (Associated Gene Name)')}

    validation
  end

  dep :validation_dataset
  input :only_consensus, :boolean, "Use only consensus sentences", true
  task :aug_validation_dataset => :tsv do |only_consensus|
    validation_dataset = step(:validation_dataset).load
    validation_consensus = Rbbt.data["Astrid Validation"]["Consensus_280717_extended_data_v3.xlsx - consensus.tsv"].tsv :type => :list
    validation_HC = Rbbt.data["Astrid Validation"]["validation_384_high-confidence_ExTRI_July_2017.xlsx - validation_July_2017.tsv"].tsv :type => :list

    validation = TSV.setup({}, :key_field => "PMID:Sentence ID:TF:TG", :fields => ["Valid"], :type => :single)

    validation_dataset.through do |k,values|
      validation[k] = values["Valid"] == "Valid" ? "Valid" : "Not Valid"
    end unless only_consensus

    validation_HC.through do |k,values|
      validation[k] = values["validation"] == "V" ? "Valid" : "Not Valid"
    end unless only_consensus

    validation_consensus.through do |k,values|
      consensus = values["Consensus"]
      next if consensus == ""
      validation[k] = consensus == "V" ? "Valid" : "Not Valid"
    end

    validation
  end

  #{{{ Prediction model

  dep :aug_validation_dataset
  dep :ExTRI_counts
  dep :ExTRI_postprocess
  input :post_process, :boolean, "Filter training sentences by postprocessing rules", false
  input :test_set, :array, "Separate entries for testing", []
  task :prediction => :tsv do |post_process,test_set|
    full = step(:ExTRI_counts).load
    validation = step(:aug_validation_dataset).load

    if post_process
      post = step(:ExTRI_postprocess).load
      validation = validation.select(validation.keys & post.keys)
    end

    if test_set and test_set.any?
      validation = validation.select(validation.keys - test_set)
    end

    all_validation_keys = validation.keys
    valid_keys = validation.select("Valid" => "Valid").keys

    require 'rbbt/util/R'
    train = full.select(all_validation_keys).attach validation, :fields => ["Valid"]
    Open.write(file('train'), train.to_s)
    data = nil
    TmpFile.with_file do |file|
      train.R <<-EOF
      rbbt.require('randomForest')
      names(data) <- make.names(names(data))
      data$Valid <- as.factor(data$Valid)
      m = randomForest(Valid ~ Interaction.score + PMID.counts + Sentence.counts + Sentence.pairs + Sentence.pair.density, data=data)
      #m = randomForest(Valid ~ Interaction.score + PMID.counts + Sentence.counts, data=data)
      save(m, file='#{file}')
      EOF

      data = full.R <<-EOF
      rbbt.require('randomForest')
      names(data) <- make.names(names(data))
      load('#{file}')
      predictions = predict(m, data)
      predictions = predictions[!is.na(predictions)]
      data = data.frame(predictions)
      rownames(data) <- names(predictions)
      EOF
    end

    data.select("predictions" => "Valid")
  end

  #{{{ Add confidence to dataset
  
  dep :ExTRI_counts, :compute => :produce
  dep :prediction, :compute => :produce
  dep :threshold, :compute => :produce
  desc <<-EOF
Takes the ExTRI_counts file and adds the prediction and threshold confidence. 

Both confidence calls are force to Low if the target gene is a signal transduction element.
  EOF
  task :ExTRI_confidence => :tsv  do
    predicted = Set.new step(:prediction).load.keys
    thresholded = Set.new step(:threshold).load.keys
    signal_transd = Rbbt.data["signal_transd.list"].list

    tsv = step(:ExTRI_counts).load
    tsv.with_monitor do
    tsv.add_field "Prediction confidence" do |k,v|
      tg = v["Target Gene (Associated Gene Name)"]
      if signal_transd.include? tg
        "Low"
      else
        predicted.include?(k) ? "High" : "Low"
      end
    end

    score, pmids, sentences  = step(:threshold).inputs.values_at :score, :pmids, :sentences
    criteria = "score: #{score}; min # PMID: #{pmids}; min # Sentences in PMID: #{sentences}"
    tsv.add_field "Threshold confidence (#{criteria})" do |k,v|
      tg = v["Target Gene (Associated Gene Name)"]
      if signal_transd.include? tg
        "Low"
      else
        thresholded.include?(k) ? "High" : "Low"
      end
    end
    end

    tsv
  end

  dep :ExTRI_confidence
  input :confidence, :select, "Confidence criteria", "Predicted", :select_options => ["Predicted", "Threshold"]
  desc <<-EOF
Assigns confidence for every ExTRI triplet (TF:TG:PMID) based on the best confidence call for it (best sentence confidence)
  EOF
  task :triplet_confidence => :tsv do |confidence|
    tsv = step(:ExTRI_confidence).load
    confidence_field = confidence == "Predicted" ? "Prediction confidence" : tsv.fields.select{|f| f =~ /Threshold/}.first

    res =  TSV.setup({}, :key_field => "TF:TG:PMID", :fields => [confidence_field], :type => :single)
    tsv.column(confidence_field).through do |k,confidence|
      pmid,sent,tf,tg = k.split(":")
      triplet = [tf,tg,pmid] * ":"
      res[triplet] = confidence if res[triplet].nil? or res[triplet] == "Low"
    end
    res
  end

  #dep :validation_dataset
  #dep :sentence_coverage_full_NER_counts
  #task :predicted_pairs => :array do 
  #  full = step(:sentence_coverage_full_NER_counts).load
  #  validation = step(:validation_dataset).load

  #  all_validation_keys = validation.keys
  #  valid_keys = validation.select("Valid" => true).keys

  #  require 'rbbt/util/R'
  #  train = full.select(all_validation_keys).attach validation, :fields => ["Valid"]
  #  Open.write(file('train'), train.to_s)
  #  data = nil
  #  TmpFile.with_file do |file|
  #    train.R <<-EOF
  #    library(randomForest)
  #    names(data) <- make.names(names(data))
  #    data$Valid <- as.factor(data$Valid)
  #    m = randomForest(Valid ~ Interaction.score + PMID.counts + Sentence.counts, data=data)
  #    save(m, file='#{file}')
  #    EOF

  #    data = full.R <<-EOF
  #    library(randomForest)
  #    names(data) <- make.names(names(data))
  #    load('#{file}')
  #    predictions = predict(m, data)
  #    predictions = predictions[!is.na(predictions)]
  #    data = data.frame(predictions)
  #    rownames(data) <- names(predictions)
  #    EOF
  #  end

  #  data.select("predictions" => "Valid").keys.collect{|k| k.split(":").values_at(2,3)*":"}
  #end

  #dep :predicted_pairs
  #dep :threshold_pairs, :pmids => 1000 do |jobname,options|
  #  pmids = [2, 3, 1000]
  #  o = options.dup
  #  pmids.collect do |pmid|
  #    o[:pmids] = pmid
  #    ExTRI.job(:threshold_pairs, jobname, o)
  #  end
  #end
  #task :pair_analysis => :tsv do
  #  pred = step(:predicted_pairs).load
  #  tsv = TSV.setup({}, :key_fields => "TF:TG", :fields => [], :type => :list)

  #  pred_tsv = TSV.setup(pred, :key_fields => "TF:TG", :fields => [], :type => :list)
  #  pred_tsv.add_field "Prediction" do 
  #    "Prediction"
  #  end
  #  tsv.attach pred_tsv, :complete => true

  #  dependencies.each do |job|
  #    next if job.path.include? 'predicted_pairs'
  #    pmids = job.recursive_inputs[:pmids]
  #    t = TSV.setup(job.load, :key_fields => "TF:TG", :fields => [], :type => :list)
  #    t.add_field "Threshold PIMD > #{pmids}" do 
  #      "Threshold #{pmids}"
  #    end
  #    tsv.attach t, :complete => true
  #  end

  #  tsv
  #end
end
