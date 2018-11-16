module ExTRI

  input :salvage, :boolean, "Salvage some non-TFClass TFs", false
  task :flagged_tfs_TFCheckpoint => :array do |salvage|

    cp = TFCheckpoint.tfs.tsv

    new_cp = cp.annotate({})

    translation = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :order => true, :persist => true
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

    fnl_fs = Open.open(ExTRI.TF_full) do |io| Misc.fixutf8(CMD.cmd("cut -f 2 | sort -u", :in => io).read).split("\n") end
    fnl_fs - good_tfs - ["NFKB", "AP1"]
  end

  input :salvage, :boolean, "Salvage some non-TFClass TFs", false
  task :flagged_tfs_GTRD => :array do |salvage|
    require 'rbbt/sources/GTRD'

    ens = GTRD.tfClass.list

    translation = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :order => true, :persist => true

    good_tfs = translation.values_at(*ens).compact

    fnl_fs = Open.open(ExTRI.TF_full) do |io| Misc.fixutf8(CMD.cmd("cut -f 2 | sort -u", :in => io).read).split("\n") end
    fnl_fs - good_tfs - ["NFKB", "AP1"]
  end

  dep :flagged_tfs_GTRD
  task :flagged_tfs => :array do
    dependencies.first.load
  end

  task :ExTRI_clean  => :tsv do
    p = TSV::Parser.new ExTRI.TF_full.find

    translation = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :order => true, :persist => true
    dumper = TSV::Dumper.new :key_field => "PMID:Sentence ID:TF:TG", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "Interaction score", "Sentence"], :type => :list, :namespace => ExTRI.organism
    dumper.init
    tf_pos = 0
    TSV.traverse ExTRI.TF_full.find, :type => :list, :into => dumper do |k,v|
      v = v.values_at 0,1,2,3

      if ! %w(AP1 NFKB).include?(v[0])
        tf = translation[v[0]] || v[0]
        Log.low "Translate #{v[0]} => #{ tf }" if v[0] != tf
        v[0] = tf 
      end

      if ! %w(AP1 NFKB).include?(v[1])
        tg = translation[v[1]] || v[1]
        Log.low "Translate #{v[1]} => #{ tg }" if v[1] != tg
        v[1] = tg 
      end

      key = [k, v[0], v[1]] * ":"
      [key,v]
    end

    Misc.sort_stream(dumper.stream)
  end


  helper :sentence_contains do |sentence,pattern|
    parts = pattern.split(/[&!]/)
    matches = parts.collect do |part|
      sent = part.downcase == part ? sentence.downcase : sentence
      if part =~ /^\w+$/
        sent.split(/[^\w]/).include? part
      else
        sent =~ /\W#{part}\W/i
      end
    end
    if pattern.include? "&"
      matches.inject(true){|acc,e| acc = acc && e}
    elsif pattern.include? "|"
      matches.inject(false){|acc,e| acc = acc || e}
    else
      matches.first
    end
  end

  dep :ExTRI_clean
  task :ExTRI_postprocess => :tsv do
    tsv = step(:ExTRI_clean).load

    TSV.traverse Rbbt.root.data["post_process_rules.tsv"].find(:lib), :type => :array do |line|
      next if line =~ /^#/
      rtf,rtg,has,hasnot,exact = line.split(",").collect{|v| v.empty? ? nil : v}

      gene = [rtf,rtg].compact.first

      rejects = []

      tsv = tsv.select do |k,v|
        tf, tg, *rest = v
        sentence = v["Sentence"]

        reject = false
        if (rtf and tf =~ /^#{rtf}$/) or (rtg and tg =~ /^#{rtg}$/)
          if has 
            reject = true if sentence_contains(sentence, has)
          else
            reject = true
          end

          if hasnot
            reject = false if sentence_contains(sentence, hasnot)
          end

          if exact == 'true'
            reject = false if sentence_contains(sentence, gene)
          end
        end

        rejects << [k, v].flatten * "\t" if reject
        ! reject
      end

      Open.write(file(gene), rejects * "\n")
      log gene, "Removed #{gene}. Has: '#{has}'; has not: '#{hasnot}'; rescue exact HGNC: '#{exact}' - #{rejects.length} removed"
    end

    log :HDACS_list, "Remove all HDACs" do
      rejects = []
      hdacs = Rbbt.root.data["hdacs.list"].list

      tsv = tsv.select do |k,v|
        tf, tg, *rest = v

        reject = hdacs.include? tg

        rejects << [k, v].flatten * "\t" if reject
        ! reject
      end

      Open.write(file("HDACS_list"), rejects * "\n")
    end

    tsv
  end

  dep :ExTRI_postprocess
  input :skip_post_process, :boolean, "Count entities including sentences removed by post_processing", false
  task :ExTRI_counts => :tsv do |skip_post_process|

    pmid_counts = TSV.setup({}, :key_field => 'Pair', :fields => ["Counts"], :type => :single)
    sentence_counts = TSV.setup({}, :key_field => 'Triplet', :fields => ["Counts"], :type => :single)
    sentence_pairs = TSV.setup({}, :key_field => 'Sentence', :fields => ["Counts"], :type => :single)

    dep = skip_post_process ? step(:ExTRI_postprocess).step(:ExTRI_clean) : step(:ExTRI_postprocess)

    TSV.traverse dep do |k,values|
      tf, tg = values.values_at(0,1)
      sentence = k.split(":").values_at(0,1) * ":"
      pair = [tf, tg] * ":"
      triplet = [tf, tg, k.split(":").first] * ":"
      pmid_counts[pair] ||= 0
      pmid_counts[pair] += 1
      sentence_counts[triplet] ||= 0
      sentence_counts[triplet] += 1
      sentence_pairs[sentence] ||= 0
      sentence_pairs[sentence] += 1
    end

    dumper = TSV::Dumper.new :key_field => "PMID:Sentence ID:TF:TG", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "Interaction score", "Sentence", "PMID counts", "Sentence counts", "Sentence pairs", "Sentence length", "Sentence pair density"], :type => :list, :namespace => ExTRI.organism
    dumper.init
    TSV.traverse dep, :into => dumper do |k,values|
      tf, tg = values.values_at(0,1)
      pair = [tf, tg] * ":"
      triplet = [tf, tg, k.split(":").first] * ":"
      sentence = k.split(":").values_at(0,1) * ":"
      pmid_c = pmid_counts[pair] 
      sentence_c = sentence_counts[triplet]
      sentence_p = sentence_pairs[sentence]
      [k, values + [pmid_c, sentence_c, sentence_p, sentence.length, sentence_p.to_f / sentence.length]]
    end
  end

  #dep :ExTRI_postprocess
  #task :ExTRI_upstream_regulators => :tsv do
  #  upstream = Rbbt.data["hdacs_etc.list"].find(:lib).list
  #  TSV.traverse step(:ExTRI_postprocess), :type => :array, :into => :stream do |line|
  #    tg = line.split("\t")[2] 
  #    next unless upstream.include? tg
  #    line
  #  end
  #end

  #dep :flagged_tfs
  #task :sentence_coverage_subset => :tsv do

  #  flagged_tfs = step(:flagged_tfs).load

  #  id_file = Organism.identifiers(ExTRI.organism)

  #  encode = ExTRI.Encode.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
  #  goa = ExTRI.GOA.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
  #  intact = ExTRI.Intact.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip

  #  tfacts = TFacts.tf_tg.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :merge => true, :zipped => true).unzip
  #  trrust = TRRUST.tf_tg.tsv(:merge => true).unzip
  #  htri = HTRI.tf_tg.tsv(:merge => true).unzip

  #  flagged = ExTRI.TFacts_flagged_articles.list
  #  tfacts.add_field "Confidence" do |tf,values|
  #    sign,species,source,pmids = values
  #    (source.downcase == "pubmed" and (pmids.split(";") - flagged).empty?) ? "Low" : "High"
  #  end

  #  tsv = ExTRI.TF.tsv(:merge => true).select("TF Associated Gene Name").unzip(0, true, ":", false).unzip("TF Associated Gene Name", true, ":", false).unzip("TG Associated Gene Name", true, ":", false).to_list

  #  tsv = attach_db tsv, htri, "HTRI"
  #  tsv = attach_db tsv, trrust, "TRRUST"
  #  tsv = attach_db tsv, tfacts, "TFacts"
  #  tsv = attach_db tsv, encode, "Encode"
  #  tsv = attach_db tsv, goa, "GOA"
  #  tsv = attach_db tsv, intact, "Intact"

  #  tsv
  #end
end
