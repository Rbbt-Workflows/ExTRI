module FNL

  input :salvage, :boolean, "Salvage some non-TFClass TFs", false
  task :flagged_tfs_TFCheckpoint => :array do |salvage|

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

    fnl_fs = Open.open(FNL.TF_full) do |io| Misc.fixutf8(CMD.cmd("cut -f 2 | sort -u", :in => io).read).split("\n") end
    fnl_fs - good_tfs - ["NFKB", "AP1"]
  end

  input :salvage, :boolean, "Salvage some non-TFClass TFs", false
  task :flagged_tfs_GDRD => :array do |salvage|
    require 'rbbt/sources/GTRD'

    ens = GTRD.tfClass.list

    translation = Organism.identifiers(FNL.organism).index :target => "Associated Gene Name", :order => true, :persist => true

    good_tfs = translation.values_at(*ens).compact

    fnl_fs = Open.open(FNL.TF_full) do |io| Misc.fixutf8(CMD.cmd("cut -f 2 | sort -u", :in => io).read).split("\n") end
    fnl_fs - good_tfs - ["NFKB", "AP1"]
  end

  dep :flagged_tfs_GDRD
  task :flagged_tfs => :array do
    dependencies.first.load
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
  task :FNL_postprocess => :tsv do
    tsv = step(:FNL_clean).load

    log :LPS, "Remove LPS" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"]
        select = if tf == "IRF6" or tg == "IRF6"
                   not sentence.include?("LPS") 
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:LPS), rejects)
    end

    log :GNAS, "Remove GNAS" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tf == "GNAS" or tg == "GNAS"
                   not (sentence.include?("promoter") and sentence.include?("p2"))
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:GNAS), rejects)
    end

    log :ADAM, "Remove ADAM" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tf == "ADAM2" or tg == "ADAM2"
                   sentence.include?("adam")
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:GNAS), rejects)
    end

    log :FOXC1, "Remove FOXC1" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tf == "FOXC1" or tg == "FOXC1"
                   not sentence.include?("chip")
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:FOXC1), rejects)
    end


    log :CAT, "Remove CAT" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tf == "CAT" or tg == "CAT"
                   sentence.include?("catalase")
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:CAT), rejects)
    end

    log :HDAC, "Remove HDAC" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg.include?("HDAC")
                   false
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:HDAC), rejects)
    end

    log :ESR1, "Remove ESR1" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tf == "ESR1" or tg == "ESR1"
                   not sentence.include?("stress")
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:ESR1), rejects)
    end

    log :MAPK, "Remove MAPK" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg =~ /^MAPK\d*$/
                   false
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:MAPK), rejects)
    end

    log :PKC, "Remove PKC" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg =~ /^PKC\d*$/
                   false
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:PKC), rejects)
    end

    log :AKT, "Remove AKT" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg =~ /^AKT\d*$/
                   false
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:AKT), rejects)
    end

    log :PI3K, "Remove PI3K" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg =~ /^PI3K\d*$/
                   false
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:PI3K), rejects)
    end


    log :CEL, "Remove CEL" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg == "CEL" or tf == "CEL"
                   not sentence.include?("cell")
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:CEL), rejects)
    end

    log :ATP11C, "Remove ATP11C" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg == "ATP11C"
                   not sentence.split(/[^\w]/).include? "ig"
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:ATP11C), rejects)
    end

    log :SUPT7L, "Remove SUPT7L" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg == "SUPT7L"
                   not sentence.split(/[^\w]/).include? "gamma"
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:SUPT7L), rejects)
    end

    log :TPM1, "Remove TPM1" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg == "TPM1"
                   not sentence.split(/[^\w]/).include? "alpha"
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:TPM1), rejects)
    end

    log :DLST, "Remove DLST" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg == "DLST"
                   not sentence.split(/[^\w]/).include? "e2"
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:DLST), rejects)
    end

    log :TSC1, "Remove TSC1" do

      rejects = ""
      tsv = tsv.select do |k,v| 
        tf, tg, *rest = v
        sentence = v["Sentence"].downcase
        select = if tg == "TSC1"
                   not sentence.split(/[^\w]/).include? "suppressor"
                 else
                   true
                 end
        rejects << [k, v].flatten * "\t" << "\n" if not select
        select
      end

      Open.write(file(:TSC1), rejects)
    end

    tsv
  end

  dep :FNL_postprocess
  task :FNL_counts => :tsv do

    pmid_counts = TSV.setup({}, :key_field => 'Pair', :fields => ["Counts"], :type => :single)
    sentence_counts = TSV.setup({}, :key_field => 'Triplet', :fields => ["Counts"], :type => :single)

    TSV.traverse step(:FNL_postprocess) do |k,values|
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
    TSV.traverse step(:FNL_postprocess), :into => dumper do |k,values|
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
