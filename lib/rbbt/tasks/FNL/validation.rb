module FNL

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
  dep :FNL_clean
  task :greco_format => :text do 
    tsv = step(:validation_dataset).load.select("Valid" => "Valid")
    fixed = step(:FNL_clean).load

    name2ens = Organism.identifiers("Hsa/feb2014").index :persist => true
    hashes = []
    Log.tsv fixed
    tsv.through do |pair,values|
      iii [pair, fixed.include?(pair)]
      sentence = fixed[pair]["Sentence"]
      pmid, n, tf, tg = pair.split(":")

      ens_tf = name2ens[tf]
      ens_tg = name2ens[tg]

      url_tf = "http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=#{ens_tf}"
      url_tg = "http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=#{ens_tg}"
      
      
      hashes << {
        "ext_id": pmid, #pmid
        "provider": "FNL", 
        "anns": [
          {
            "exact": sentence,
            "section": "abstract",
            "tags": [
              {"name": tf,
               "uri": url_tf },
               {"name": tg,
                "uri": url_tg }
            ]
          }
        ]
      }
    end


    hashes.to_json
  end

  dep :validation_dataset
  dep :sentence_coverage_NER
  input :target, :integer, "Target number of sentences", 2000
  task :extended_validation_dataset => :tsv do |target|
    dataset = step(:validation_dataset).load
    tsv = step(:sentence_coverage_NER).load
    current = dataset.keys
    missing = target - current.length
    if missing > 0 
      dataset_fields = dataset.fields.collect{|f| f.gsub(/[()]/,'').sub("Transcription Factor", "TF").sub("Target Gene", "TG") }
      tsv_fields = tsv.fields.collect{|f| f.gsub(/[()]/,'').sub("Transcription Factor", "TF").sub("Target Gene", "TG") }

      new_ids = (tsv.keys - current).shuffle[0..missing-1]
      new_ids.each do |id|
        values = tsv[id]
        dataset_values = [nil] * dataset_fields.length
        tsv_fields.zip(values).each do |field,value|
          pos = dataset_fields.index field.gsub(/[()]/,'')
          dataset_values[pos] = value if pos
        end

        dataset_values[dataset_fields.index("Set")] = "Extended"

        dataset[id] = dataset_values
      end
      dataset
    else
      dataset
    end

    keys = dataset.keys
    new_keys = keys[0..99]
    new_keys.concat keys[100..-1].shuffle

    new = dataset.annotate({})
    new_keys.each do |k|
      new[k] = dataset[k]
    end

    new
  end

  input :score, :float, "Score threshold", 1.6
  input :pmids, :integer, "Min number of PMIDS", 1000
  input :sentences, :integer, "Min number of sentences", 2
  dep :FNL_counts, :compute => :produce
  task :threshold => :tsv do |score, pmids, sentences|
    tsv = step(:FNL_counts).load
    selected = []
    selected.concat tsv.select("Interaction score"){|s| s.to_f >= score }.keys
    selected.concat tsv.select("PMID counts"){|s| s.to_f >= pmids }.keys
    selected.concat tsv.select("Sentence counts"){|s| s.to_f >= sentences }.keys

    tsv.select(selected)
  end

  dep :threshold 
  dep :validation_dataset
  task :threshold_evaluation => :yaml do 
    full = step(:threshold).step(:FNL_counts).load
    tsv = step(:threshold).load

    validation = step(:validation_dataset).load

    all_validation_keys = validation.keys

    require 'rbbt/util/R'
    train = full.select(all_validation_keys)
    train = train.attach validation, :fields => ["Valid"]
    data = nil
    TmpFile.with_file do |file|
      train.R <<-EOF
      library(randomForest)
      names(data) <- make.names(names(data))
      data$Valid <- as.factor(data$Valid)
      m = randomForest(Valid ~ Interaction.score + PMID.counts + Sentence.counts, data=data)
      save(m, file='#{file}')
      EOF

      data = tsv.R <<-EOF
      library(randomForest)
      names(data) <- make.names(names(data))
      load('#{file}')
      predictions = predict(m, data)
      data = data.frame(predictions)
      EOF
    end

    res = Misc.counts(data.column("predictions").values.flatten)
    res["percent"] = 100 * res["Valid"] / (res["Valid"] + res["Not Valid"])
    res
  end

  dep :threshold
  task :threshold_pairs => :array do 
    step(:threshold).load.keys.collect{|k| k.split(":").values_at(2,3) * ":" }.uniq
  end

  dep :validation_dataset
  dep :FNL_counts
  task :prediction => :tsv do 
    full = step(:FNL_counts).load
    validation = step(:validation_dataset).load

    all_validation_keys = validation.keys
    valid_keys = validation.select("Valid" => true).keys

    require 'rbbt/util/R'
    train = full.select(all_validation_keys).attach validation, :fields => ["Valid"]
    Open.write(file('train'), train.to_s)
    data = nil
    TmpFile.with_file do |file|
      train.R <<-EOF
      library(randomForest)
      names(data) <- make.names(names(data))
      data$Valid <- as.factor(data$Valid)
      m = randomForest(Valid ~ Interaction.score + PMID.counts + Sentence.counts, data=data)
      save(m, file='#{file}')
      EOF

      data = full.R <<-EOF
      library(randomForest)
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

  dep :FNL_counts, :compute => :produce
  dep :prediction, :compute => :produce
  dep :threshold, :compute => :produce
  task :FNL_confidence => :tsv  do
    predicted = Set.new step(:prediction).load.keys
    thresholded = Set.new step(:threshold).load.keys

    tsv = step(:FNL_counts).load
    tsv.with_monitor do
    tsv.add_field "Prediction confidence" do |k,v|
      predicted.include?(k) ? "High" : "Low"
    end

    score, pmids, sentences  = step(:threshold).inputs.values_at :score, :pmids, :sentences
    criteria = "score: #{score}; min # PMID: #{pmids}; min # Sentences in PMID: #{sentences}"
    tsv.add_field "Threshold confidence (#{criteria})" do |k,v|
      thresholded.include?(k) ? "High" : "Low"
    end
    end

    tsv
  end

  dep :FNL_confidence
  input :confidence, :select, "Confidence criteria", "Predicted", :select_options => ["Predicted", "Threshold"]
  task :triplet_confidence => :tsv do |confidence|
    tsv = step(:FNL_confidence).load
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
  #    FNL.job(:threshold_pairs, jobname, o)
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
