module FNL

  dep :threshold 
  dep :aug_validation_dataset
  task :threshold_evaluation => :yaml do 
    full = step(:threshold).step(:FNL_counts).load
    tsv = step(:threshold).load
    validation = step(:aug_validation_dataset).load

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

  dep :aug_validation_dataset
  dep :FNL_confidence
  input :test_set, :array, "Separate entries for testing", []
  task :FNL_confidence_evaluation => :tsv do |test_set|
    vali = step(:aug_validation_dataset).load
    conf = step(:FNL_confidence).load

    if test_set
      vali = vali.select(test_set)
    end

    valid = vali.select("Valid" => "Valid").keys
    novalid = vali.keys - valid
    
    all = vali.keys

    cvalid = conf.select("Prediction confidence" => "High").keys & all
    cnovalid = (conf.keys - cvalid) & all

    tp = cvalid & valid
    fp = cvalid & novalid
    tn = cnovalid & novalid
    fn = cnovalid & valid

    tsv = TSV.setup({}, :key_field => "Statistic", :fields => ["Value"], :cast => :to_f, :type => :single)
    tsv["TP"] = tp.length
    tsv["FP"] = fp.length
    tsv["TN"] = tn.length
    tsv["FN"] = fn.length
    tsv["TPR"] = tp.length.to_f / cvalid.length
    tsv["TNR"] = tn.length.to_f / cnovalid.length


    tsv
    
  end

  dep :aug_validation_dataset, :compute => :produce
  input :cv_times, :integer, "Number of CV folds", 5
  dep :FNL_confidence_evaluation, :test_set => :placeholder, :compute => :bootstrap do |jobname,options,dependencies|
    keys = dependencies.first.run.keys

    cv_times = options[:cv_times]
    size = keys.length / cv_times
    sets = []
    cv_times.times.each do 
      sets << keys.shuffle[0..size-1]
    end

    sets.collect do |set|
      FNL.job(:FNL_confidence_evaluation, jobname, options.merge({:test_set => set}))
    end
  end
  task :FNL_confidence_CV => :tsv do
    res = nil
    i = 0
    dependencies.collect do |dep|
      next unless dep.task_name.to_s == "FNL_confidence_evaluation"
      i += 1
      tsv = dep.load
      tsv.fields = tsv.fields.collect{|f| f + " (#{i})" }
      if res.nil? 
        res = tsv
      else
        res = res.attach tsv
      end
    end
    res.add_field "Average" do |k,values|
      Misc.mean(values)
    end
    res
  end

  dep :FNL_confidence, :only_consensus => :placeholder do |jobname, options|
    jobs = []
    jobs << FNL.job(:FNL_confidence, jobname, options.merge(:only_consensus => true))
    jobs << FNL.job(:FNL_confidence, jobname, options.merge(:only_consensus => false))
    jobs
  end
  task :FNL_confidence_differences => :tsv do
    only = dependencies.first.load
    all = dependencies.last.load

    only_valid = only.select("Prediction confidence" => "High").keys 
    all_valid = all.select("Prediction confidence" => "High").keys 
    
    tsv = TSV.setup({}, :key_field => "Statistic", :fields => ["value"], :type => :list, :cast => :to_i)
    tsv["Common"] = (only_valid & all_valid).length
    tsv["Only consensus"] = (only_valid - all_valid).length
    tsv["All"] = (all_valid - only_valid).length
    tsv 
  end



end
