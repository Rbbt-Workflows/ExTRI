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

    validation = ExTRI.update_tsv_symbols(validation)

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

    validation = ExTRI.update_tsv_symbols(validation)

    validation
  end

  dep :aug_validation_dataset
  task :liv_validation => :tsv do
    validation_dataset = step(:aug_validation_dataset).load
    validation_liv = TSV.excel(Rbbt.data["Astrid and Liv validation"]["curated_sentences_compiled_270121.xlsx"].find, :type => :list)

    validation_liv.through do |k,values|
      val = values["evaluation"] == 'v'
      validation_dataset[k] = val ? "Valid" : "Not Valid"
    end

    validation_dataset = ExTRI.update_tsv_symbols(validation_dataset)

    validation_dataset
  end

  task :NTNU_curation_validation => :tsv do
    tsv = Rbbt.data.NTNU_Curation.current_stack.tsv :fields => ["Valid"]
    tsv.key_field = "PMID:Sentence ID:TF:TG"
    tsv.process "Valid" do |v|
      v.to_s == 'true' ? "Valid" : "Not Valid"
    end
    tsv
  end

  dep :NTNU_curation_validation
  task :validation_dataset => :tsv do 
    tsv = dependencies.first.load
    tsv.key_field = "PMID:Sentence ID:TF:TG"
    ExTRI.update_tsv_symbols tsv
  end

  #{{{ Prediction model

  dep :validation_dataset
  dep :ExTRI_counts
  dep :ExTRI_postprocess
  input :post_process, :boolean, "Filter training sentences by postprocessing rules", false
  input :test_set, :array, "Separate entries for testing", []
  task :prediction => :tsv do |post_process,test_set|
    full = step(:ExTRI_counts).load
    validation = step(:validation_dataset).load

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
    train = full.select(all_validation_keys)
    train.attach validation, :fields => ["Valid"]
    Open.write(file('train'), train.to_s)
    data, predictions, scores = nil
    TmpFile.with_file do |file|
      train.R <<-EOF
      rbbt.require('randomForest')
      names(data) <- make.names(names(data))
      data$Valid <- as.factor(data$Valid)
      m = randomForest(Valid ~ Interaction.score + PMID.counts + Sentence.counts + Sentence.pairs + Sentence.pair.density, data=data)
      #m = randomForest(Valid ~ Interaction.score + PMID.counts + Sentence.counts, data=data)
      save(m, file='#{file}')
      EOF

      predictions = full.R <<-EOF
      rbbt.require('randomForest')
      names(data) <- make.names(names(data))
      load('#{file}')
      predictions = predict(m, data)
      predictions = predictions[!is.na(predictions)]
      data = data.frame(predictions)
      rownames(data) <- names(predictions)
      EOF

      scores = full.R <<-EOF
      rbbt.require('randomForest')
      names(data) <- make.names(names(data))
      load('#{file}')
      predictions = predict(m, data)
      probs = predict(m, data, type='prob')
      probs = probs[!is.na(predictions),]
      predictions = predictions[!is.na(predictions)]
      data = data.frame(probs)
      rownames(data) <- names(predictions)
      EOF
    end
    Log.tsv predictions
    Log.tsv scores

    data = predictions.to_list.attach(scores, :fields => ["Valid"])
    data.fields = %w(Prediction Probability)

    data
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
    prediction_scores = step(:prediction).load.slice(["Probability"])
    predicted = Set.new step(:prediction).load.select("Prediction" => "Valid").keys
    thresholded = Set.new step(:threshold).load.keys
    signal_transd = Rbbt.data["signal_transd.list"].list

    tsv = step(:ExTRI_counts).load

    tsv.add_field "Prediction confidence" do |k,v|
      tg = v["Target Gene (Associated Gene Name)"]
      if signal_transd.include? tg
        "Low"
      else
        predicted.include?(k) ? "High" : "Low"
      end
    end

    tsv.add_field "Prediction confidence (score)" do |k,v|
      prediction_scores[k]
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

    tsv.add_field "Auto-regulation" do |pair,values|
      (values[0] == values[1]) ? "Auto-regulation" : ""
    end

    tsv.process "Prediction confidence" do |v,k,l|
      if l["Auto-regulation"] == "Auto-regulation"
        "Low"
      else
        v
      end
    end

    tsv
  end

  dep :ExTRI_confidence
  input :confidence, :select, "Confidence criteria", "Prediction", :select_options => ["Prediction", "Threshold"]
  desc <<-EOF
Assigns confidence for every ExTRI triplet (TF:TG:PMID) based on the best confidence call for it (best sentence confidence)
  EOF
  task :triplet_confidence => :tsv do |confidence|
    tsv = step(:ExTRI_confidence).load
    confidence_field = confidence == "Prediction" ? "Prediction confidence" : tsv.fields.select{|f| f =~ /Threshold/}.first

    res =  TSV.setup({}, :key_field => "TF:TG:PMID", :fields => [confidence_field], :type => :single)
    tsv.column(confidence_field).through do |k,confidence|
      pmid,sent,tf,tg = k.split(":")
      triplet = [tf,tg,pmid] * ":"
      res[triplet] = confidence if res[triplet].nil? or res[triplet] == "Low"
    end
    res
  end

  dep :ExTRI_confidence
  task :ExTRI_final => :tsv do
    invalid = Rbbt.root.data["NTNU_invalid"].list
    step(:ExTRI_confidence).load.select(invalid, true)
  end


end
