module FNL

  input :salvage, :boolean, "Salvage some non-TFClass TFs", false
  task :flagged_tfs => :array do |salvage|

    cp = TFCheckpoint.tfs.tsv

    new_cp = cp.annotate({})

    translation = Organism.identifiers(FNL.organism).index :target => "Associated Gene Name", :order => true, :persist => true
    TSV.traverse cp, :into => new_cp do |tf, values|
      new_tf = translation[tf] 
      if new_tf.nil?
        Log.error "Gene name not found from TFCheckpoint: " + tf
        new_tf = tf
      end
      Log.info "TFCheckpoint Translation " << [tf, new_tf] * " => " if tf != new_tf
      [new_tf, values]
    end

    cp = new_cp

    tf_class_human = cp.select("TFClass_human"){|v| v.include?("TFclass")}.keys.flatten.uniq.compact
    tf_class_mouse = cp.select("TFclass_mouse"){|v| v.include?("TFclass_mouse")}.keys.flatten.uniq.compact

    tf_class = (tf_class_human + tf_class_mouse).uniq

    if salvage
      good_tfs = (tf_class + Rbbt.data.TFclass_extended.salvage.list).uniq
    else
      good_tfs = tf_class
    end

    fnl_fs = Misc.fixutf8(CMD.cmd("cut -f 2 '#{FNL.TF_full.find}' | sort -u").read).split("\n")
    fnl_fs - good_tfs - ["NFKB", "AP1"]
  end

  dep :flagged_tfs
  task :FNL_clean  => :tsv do
    flagged_tfs = step(:flagged_tfs).load
    p = TSV::Parser.new FNL.TF_full.find

    dumper = TSV::Dumper.new :key_field => "PMID:Sentence ID:TF:TG", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "Interaction score", "Sentence"], :type => :list, :namespace => FNL.organism
    dumper.init
    tf_pos = 0
    TSV.traverse FNL.TF_full.find, :type => :list, :into => dumper do |k,v|
      v = v.values_at 0,1,2,3
      next if flagged_tfs.include? v[tf_pos]
      key = [k, v[0], v[1]] * ":"
      [key,v]
    end

    Misc.sort_stream(dumper.stream)
  end

  dep :FNL_clean
  task :FNL_counts => :tsv do

    pmid_counts = TSV.setup({}, :key_field => 'Pair', :fields => ["Counts"], :type => :single)
    sentence_counts = TSV.setup({}, :key_field => 'Triplet', :fields => ["Counts"], :type => :single)

    TSV.traverse step(:FNL_clean) do |k,values|
      tf, tg = values.values_at(0,1)
      pair = [tf, tg] * ":"
      triplet = [tf, tg, k.split(":").first] * ":"
      pmid_counts[pair] ||= 0
      pmid_counts[pair] += 1
      sentence_counts[triplet] ||= 0
      sentence_counts[triplet] += 1
    end

    dumper = TSV::Dumper.new :key_field => "PMID:Sentence ID:TF:TG", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "Interaction score", "Sentence", "PMID counts", "Sentence counts"], :type => :list, :namespace => FNL.organism
    dumper.init
    TSV.traverse step(:FNL_clean), :into => dumper do |k,values|
      tf, tg = values.values_at(0,1)
      pair = [tf, tg] * ":"
      triplet = [tf, tg, k.split(":").first] * ":"
      pmid_c = pmid_counts[pair] 
      sentence_c = sentence_counts[triplet]
      [k, values + [pmid_c, sentence_c]]
    end
  end


  #dep :flagged_tfs
  #task :sentence_coverage_subset => :tsv do

  #  flagged_tfs = step(:flagged_tfs).load

  #  id_file = Organism.identifiers(FNL.organism)

  #  encode = FNL.Encode.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
  #  goa = FNL.GOA.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
  #  intact = FNL.Intact.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip

  #  tfacts = TFacts.tf_tg.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :merge => true, :zipped => true).unzip
  #  trrust = TRRUST.tf_tg.tsv(:merge => true).unzip
  #  htri = HTRI.tf_tg.tsv(:merge => true).unzip

  #  flagged = FNL.TFacts_flagged_articles.list
  #  tfacts.add_field "Confidence" do |tf,values|
  #    sign,species,source,pmids = values
  #    (source.downcase == "pubmed" and (pmids.split(";") - flagged).empty?) ? "Low" : "High"
  #  end

  #  tsv = FNL.TF.tsv(:merge => true).select("TF Associated Gene Name").unzip(0, true, ":", false).unzip("TF Associated Gene Name", true, ":", false).unzip("TG Associated Gene Name", true, ":", false).to_list

  #  tsv = attach_db tsv, htri, "HTRI"
  #  tsv = attach_db tsv, trrust, "TRRUST"
  #  tsv = attach_db tsv, tfacts, "TFacts"
  #  tsv = attach_db tsv, encode, "Encode"
  #  tsv = attach_db tsv, goa, "GOA"
  #  tsv = attach_db tsv, intact, "Intact"

  #  tsv
  #end
end
