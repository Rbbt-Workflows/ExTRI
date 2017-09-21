module FNL

  dep :sentence_coverage_full_NER
  task :validation_sample => :tsv do
    tsv = step(:sentence_coverage_full_NER).load
    scores = tsv.column("Interaction score").values.collect{|v| v.to_f}.sort

    tsv.add_field "Rank" do |k,v|
      s = v["Interaction score"].to_f
      r = scores.index s
      r
    end


    num_bins = 10
    bin_size = (scores.length.to_f / num_bins).ceil
    bins = {}
    tsv.column("Rank").each do |k,r|
      b = (r + 1) / bin_size
      bins[b] ||= []
      bins[b] << k 
    end

    selected = []
    bins.each do |b,list|
      selected.concat list.to_a.shuffle[0..99]
    end

    selected = selected.shuffle
    new = tsv.annotate({})
    selected.each{|s| new[s] = tsv[s]}
    new.reorder(:key, ["TF Associated Gene Name", "TG Associated Gene Name", "TF Text", "TG Text", "Sentence"])
  end

  dep :sentence_coverage_full_NER
  task :validation_sample_lowscore => :tsv do
    tsv = step(:sentence_coverage_full_NER).load

    bins = {}
    tsv.column("Interaction score").each do |k,s|
      s = s.to_f
      [[0,1], [1,1.33], [1.33,1.66],[1.66,2],[2,3]].each do |range|
        start,eend = range
        next unless s > start and s <= eend
        bins[range] ||= []
        bins[range] << k
      end
    end

    selected = []
    bins.each do |b,list|
      selected.concat list.to_a.shuffle[0..99]
    end

    selected = selected.shuffle
    new = tsv.annotate({})
    selected.each{|s| new[s] = tsv[s]}
    new.reorder(:key, ["TF Associated Gene Name", "TG Associated Gene Name", "TF Text", "TG Text", "Sentence"])
  end

  dep :sentence_coverage_full_NER
  task :validation_sample_pmid => :tsv do
    tsv = step(:sentence_coverage_full_NER).load
    scores = tsv.column("Interaction score").values.collect{|v| v.to_f}.sort

    pmid_counts = TSV.setup({}, :key_field => 'Pair', :fields => ["Counts"], :type => :single)
    sentence_counts = TSV.setup({}, :key_field => 'Triplet', :fields => ["Counts"], :type => :single)

    tsv.through do |k, values|
      tf, tg = values.values_at("TF Associated Gene Name", "TG Associated Gene Name")
      pair = [tf, tg] * ":"
      triplet = k.split(":").values_at(0,2,3) * ":"
      pmid_counts[pair] ||= 0
      pmid_counts[pair] += 1
      sentence_counts[triplet] ||= 0
      sentence_counts[triplet] += 1
    end

    tsv.add_field "PMID counts" do |k,v|
      tf, tg = v.values_at("TF Associated Gene Name", "TG Associated Gene Name")
      pair = [tf, tg] * ":"
      pmid_counts[pair]
    end

    tsv.add_field "Sentence counts" do |k,v|
      triplet = k.split(":").values_at(0,2,3) * ":"
      sentence_counts[triplet]
    end

    pmid_bins = {}
    sentence_bins = {}
    tsv.column("PMID counts").each do |k,r|
      b = r
      b = -1 if b >= 7
      pmid_bins[b] ||= []
      pmid_bins[b] << k 
    end

    tsv.column("Sentence counts").each do |k,r|
      b = r
      b = -1 if b >= 3
      sentence_bins[b] ||= []
      sentence_bins[b] << k 
    end

    pmid_num = 1000 / 4
    pmid_selected = []
    pmid_bins.values_at(1,2,5,-1).each do |list|
      pmid_selected.concat list.to_a.shuffle[0..pmid_num-1]
    end

    sentence_num = 1000 / 3
    sentence_selected = []
    sentence_bins.values_at(1,2,-1).each do |list|
      sentence_selected.concat (list.to_a.shuffle - pmid_selected)[0..sentence_num-1]
    end

    selected = (pmid_selected + sentence_selected).shuffle
    new = tsv.annotate({})
    selected.each{|s| new[s] = tsv[s]}
    new.reorder(:key, ["TF Associated Gene Name", "TG Associated Gene Name", "TF Text", "TG Text", "Sentence"])
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


end
