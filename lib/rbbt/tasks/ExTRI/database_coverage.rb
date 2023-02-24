require 'rbbt/sources/TFClass'
require 'rbbt/sources/Dorothea'
module ExTRI
  AP1_SYN=%w(FOS FOSB JUN JUNB JUND FOSL1 FOSL2)
  NFKB_SYN=%w(NFKB1 NFKB2 RELA RELB)


  DATABASES=%w(ExTRI TFactS HTRI IntAct GOA TRRUST SIGNOR CytReg GEREDB)

  helper :normalize_db do |db|
    new = db.annotate({})
    translation = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :order => true, :persist => true
    
    TSV.traverse db, :into =>  new do |pair,values|
      pair = ExTRI.update_key_symbols(pair)
      values = ExTRI.update_symbols(values)
      g1, g2 = pair.split(":")
      gn1, gn2 = translation.values_at g1, g2
      #Log.error "Gene not found in translation file: " + g1 if gn1.nil?
      #Log.error "Gene not found in translation file: " + g2 if gn2.nil?
      gn1 = g1 if gn1.nil? && (NFKB_SYN + AP1_SYN).include?(g1)
      gn2 = g2 if gn2.nil? && (NFKB_SYN + AP1_SYN).include?(g2)
      next if gn1.nil? or gn2.nil?
      #Log.info "TF Translation " << [g1, gn1] * " => " if g1 != gn1
      #Log.info "TG Translation " << [g2, gn2] * " => " if g2 != gn2
      npair = [gn1, gn2] * ":"
      [npair, values]
    end

    new
  end

  helper :attach_db do |tsv, db, name|
    new = tsv.annotate({})
    new.fields += (["present"] + db.fields).collect{|f| "[#{name}] " << f }
    num_new_fields = db.fields.length

    db = normalize_db db

    TSV.traverse tsv, :into => new do |k, values|
      tf, tg = k.split(":").values_at -2, -1

      ext_tf = case tf
               when "AP1"
                 AP1_SYN
               when "NFKB"
                 NFKB_SYN
               else 
                 [tf]
               end

      ext_tg = case tg
               when "AP1"
                 AP1_SYN
               when "NFKB"
                 NFKB_SYN
               else 
                 [tg]
               end

      pairs = []
      ext_tf.each do |tf|
        ext_tg.each do |tg|
          pairs << [tf, tg] * ":"
        end
      end

      found = false
      all_values = []
      pairs.each do |pair|
        if db.include? pair
          found = true
          all_values << db[pair]
        end
      end

      if found
        extra = [name] + Misc.zip_fields(all_values).collect{|l| l*"|"} 
      else
        extra = [""] + [""] * num_new_fields
      end

      new_values = values + extra

      [k, new_values]
    end
  end

  dep :ExTRI_final, :test_set => []
  input :include_HTRI_low_conf, :boolean, "Include HTRI low confidence", false
  task :sentence_coverage_old => :tsv do |include_HTRI|
    id_file = Organism.identifiers(ExTRI.organism)

    encode = ExTRI.Encode.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip

    #goa = ExTRI.GOA.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
    goa = GO.tf_tg.tsv(:merge => true).unzip

    #intact = ExTRI.IntAct.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
    intact = IntAct.tf_tg.tsv(:merge => true).unzip

    tfacts = TFactS.tf_tg.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :merge => true, :zipped => true).unzip
    trrust = TRRUST.Hsa.tf_tg.tsv(:merge => true).unzip
    htri = HTRI.tf_tg.tsv(:merge => true).unzip(0, true)
    htri = htri.select("Confidence" => "High") unless include_HTRI
    signor = Signor.tf_tg.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => UniProt.identifiers.Hsa).unzip
    thomas = ExTRI.Thomas2015.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "sentence", "class", "details", "PMID"], :merge => true).unzip

    geredb = GEREDB.tf_tg.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "Effect", "PMID"])

    cyt_reg = CytReg.tf_cyt.tsv(:merge => true).unzip(0, true)

    flagged = ExTRI.TFactS_flagged_articles.list
    tfacts.add_field "Confidence" do |tf,values|
      sign,species,source,pmids = values
      (source.downcase == "pubmed" and (pmids.split(";") - flagged).empty?) ? "Low" : "High"
    end
  
    tsv = step(:ExTRI_final).load

    tsv = attach_db tsv, htri, "HTRI"
    tsv = attach_db tsv, trrust, "TRRUST"
    tsv = attach_db tsv, tfacts, "TFactS"
    #tsv = attach_db tsv, encode, "Encode"
    tsv = attach_db tsv, goa, "GOA"
    tsv = attach_db tsv, intact, "IntAct"
    tsv = attach_db tsv, signor, "SIGNOR"
    tsv = attach_db tsv, geredb, "GEREDB"
    tsv = attach_db tsv, cyt_reg, "CytReg"
    #tsv = attach_db tsv, thomas, "Thomas2015"

    tsv
  end

  desc <<-EOF 

List all TF:TG pairs across ExTRI and other resources along with confidence estimates and other information from those resources.

The confidence estimate for ExTRI pairs uses by default 2 PMIDs or 2 sentences or a score over 1.6.

  EOF
  dep :ExTRI_final, :pmids => 2, :sentences => 2, :score => 1.6, :test_set => []
  input :confidence, :select, "Confidence criteria", "Prediction", :select_options => ["Prediction", "Threshold"]
  input :include_HTRI_low_conf, :boolean, "Include HTRI low confidence", false
  task :pairs => :tsv do |confidence,include_HTRI|
    id_file = Organism.identifiers(ExTRI.organism)

    orig = step(:ExTRI_final).load
    signor = Signor.tf_tg.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => UniProt.identifiers.Hsa).unzip(0, true)

    pavlidis = Pavlidis.tf_tg.tsv(:merge => true)
    pavlidis = pavlidis.unzip(0, true)

    encode = ExTRI.Encode.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip

    #goa = ExTRI.GOA.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
    goa = GO.tf_tg.tsv(:merge => true).unzip

    #intact = ExTRI.IntAct.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
    intact = IntAct.tf_tg.tsv(:merge => true).unzip

    tfacts = TFactS.tf_tg.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :merge => true, :zipped => true).unzip(0, true)
    trrust = TRRUST.Hsa.tf_tg.tsv(:merge => true).unzip(0, true)

    htri = HTRI.tf_tg.tsv(:merge => true).unzip(0, true)
    htri = htri.select("Confidence" => "High") unless include_HTRI

    #thomas = ExTRI.Thomas2015.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "class", "details", "sentence", "PMID"], :merge => true).unzip
    #cp = TFCheckpoint.tfs.tsv(:merge => true)
    geredb = GEREDB.tf_tg.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "Effect", "PMID"]).unzip(0,true)

    ntnu_curated = ExTRI.NTNU_Curated.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "Sign", "PMID"], :merge => true, :type => :double).unzip(0,true)

    cyt_reg = CytReg.tf_cyt.tsv(:merge => true).unzip(0, true)

    dorotheaA = Dorothea.tf_tg.tsv(:merge => true).unzip(0, true)
    flagged = ExTRI.TFactS_flagged_articles.list
    tfacts.add_field "Confidence" do |tf,values|
      sign,species,source,pmids = values.collect{|v| v * ";"}
      (source.downcase == "pubmed" and (pmids.split(";") - flagged).empty?) ? "Low" : "High"
    end

    tsv = TSV.setup({}, :key_field => "TF:TG", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "[ExTRI] Confidence", "[ExTRI] PMID"], :type => :double, :namespace => ExTRI.organism)

    confidence_field = orig.fields.select{|f| f.include? confidence}.first
    pmids = {}
    conf = {}
    orig.through do |key,values|
      c = values[confidence_field]
      pmid, s, tf, tg = key.split(":")
      pair =  [tf,tg]
      pmids[pair] ||= []
      pmids[pair] << pmid
      conf[pair] = false if conf[pair].nil?
      conf[pair] = true unless c == "Low"
    end

    pmids.each do |pair,pmids|
      tsv[pair*":"] = [pair[0], pair[1], conf[pair] ? "High" : "Low", pmids * ";"]
    end

    tsv.add_field "[ExTRI] present" do
      "ExTRI"
    end


    [
      [htri, "HTRI"],
      [trrust, "TRRUST"],
      [tfacts, "TFactS"],
      #[encode, "Encode"], # Don't use Encode
      [goa, "GOA"],
      [intact, "IntAct"],
      [signor, "SIGNOR"],
      [cyt_reg, "CytReg"],
      [geredb, "GEREDB"],
      [ntnu_curated, "NTNU Curated"],
      [pavlidis, "Pavlidis2021"],
      [dorotheaA, "DoRothEA_A"]
      #[thomas, "Thomas2015"]
    ].each do |db,name|
      log :adding_db, name

      db.key_field = tsv.key_field

      tsv = attach_db tsv, db, name

      db = normalize_db db

      db.through do |k, values|
        next if tsv.include? k
        begin
          new_values = k.split(":") + ([""] * (tsv.fields.length - 3 - db.fields.length)) + [name] + values
          tsv[k] = new_values
        rescue
          raise $!
        end
      end
    end

    Rbbt.share.databases.ExTRI.extra_databases.glob("*.tsv").each do |file|
      name = File.basename(file).sub(/\.tsv$/,'')
      db = TSV.open(file, :merge => true).unzip

      log :adding_extra_db, name

      db.key_field = tsv.key_field

      tsv = attach_db tsv, db, name

      db = normalize_db db

      db.through do |k, values|
        next if tsv.include? k
        begin
          new_values = k.split(":") + ([""] * (tsv.fields.length - 3 - db.fields.length)) + [name] + values
          tsv[k] = new_values
        rescue
          raise $!
        end
      end
    end

    tsv.process "Transcription Factor (Associated Gene Name)" do |list,key,values|
      key.split(":").first
    end

    tsv.process "Target Gene (Associated Gene Name)" do |list,key,values|
      key.split(":").last
    end

    tsv
  end

  dep :pairs
  task :pairs_final_old => :tsv do
    tsv = step(:pairs).load
    categories = Rbbt.root.data["TFs_compExTRI_categorized_260122.tsv"].tsv :type => :single, :fields => %w(category)
    tsv.add_field "TF Category" do |pair,values|
      tf = pair.split(":").first
      category = categories[tf]
      category = "none" if category.nil?
      category.to_s.sub("_putative", "")
    end

    tfclass = TFClass.tfs.list
    tfclass << "AP1"
    tfclass << "NFKB"
    tsv.add_field "TFClass" do |pair,values|
      tf = pair.split(":").first
      (tfclass.include? tf) ? "TFClass" : ""
    end

    tsv.add_field "Auto-regulation" do |pair,values|
      (values[0] == values[1]) ? "Auto-regulation" : ""
    end

    tsv.select({"TF Category" => "none"}, true)
  end

  #dep :pairs
  #task :pairs_final_old => :tsv do
  #  tsv = step(:pairs).load
  #  tfc2 = Rbbt.root.data["TFC2_master_24012023.tsv"].tsv :type => :list, :fields => %w(lambert_2018.present TFclass_human TFclass_mouse TFclass_rat Lovering_2021.present GO:0140223.Evidence GO:0003700.Evidence GO:0003712.Evidence)
  #  tfc2["AP1"] = tfc2["JUN"]
  #  tfc2["NFKB"] = tfc2["NFKB1"]
  #  tfc2.key_field = "Transcription Factor (Associated Gene Name)"
  #  tfc2.add_field "TFClass organism" do |k,v|
  #    v.values_at(*%w(TFclass_human TFclass_mouse TFClass_rat))
  #      .compact.reject{|e| e == "NA" || e == "" }
  #      .collect{|e| e.split("_").last.split(".").first } * "+"
  #  end

  #  tfc2 = tfc2.slice(tfc2.fields - ["TFclass_mouse", "TFclass_human", "TFclass_rat"])
  #  tfc2.fields = tfc2.fields.collect{|f| f.split("_").first.split(".").first}

  #  tfc2.fields.each do |field|
  #    tfc2.process field do |v|
  #      v = field if v.include? "http"
  #      v == "NA" ? nil : v.split("_").first
  #    end
  #  end

  #  tsv.attach tfc2
  #end

  dep :pairs
  task :pairs_final => :tsv do
    tsv = step(:pairs).load
    base_path = Rbbt.share.databases.ExTRI.Feb2023_update.TF_info
    lambert_genes = base_path.Lambert_source.tsv.keys
    lovering_genes = base_path.Lovering_source_download_23022023.tsv.keys

    go_0003712_rat_uni = base_path["GO_0003712_mouse_&_rat_QuickGO_230223.tsv"].tsv(:key_field => "GENE PRODUCT ID").select("organism" => 'rat').keys
    go_0003712_mouse_uni = base_path["GO_0003712_mouse_&_rat_QuickGO_230223.tsv"].tsv(:key_field => "GENE PRODUCT ID").select("organism" => 'mouse').keys

    uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
    uni2name = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :persist => true
    gene2uniHsa = Organism.identifiers(IntAct.organism("Hsa")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniMmu = Organism.identifiers(IntAct.organism("Mmu")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniRno = Organism.identifiers(IntAct.organism("Rno")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true

    go_0003712_rat = uni2name.values_at(*uni_equivalences.values_at(*go_0003712_rat_uni).flatten.compact).compact
    go_0003712_mouse = uni2name.values_at(*uni_equivalences.values_at(*go_0003712_mouse_uni).flatten.compact).compact

    go_0003712_human = base_path["GO_0003712_human_QuickGO_230223.tsv"].tsv.keys
    go_0003700 = base_path["GO_0003700_human_QuickGO_230223.tsv"].tsv.keys
    go_0140223 = base_path["GO_0140223_human_QuickGO_230223.tsv"].tsv.keys

    %w(Lambert Lovering GO:0003700 GO:0140223).zip([lambert_genes, lovering_genes, go_0003700, go_0140223]).each do |name,genes|
      genes = ExTRI.update_symbols(genes)
      tsv.add_field name do |k,v|
        genes.include?(v.flatten.first) ? name : nil
      end
    end

    go_0003712_human = ExTRI.update_symbols(go_0003712_human)
    go_0003712_rat = ExTRI.update_symbols(go_0003712_rat)
    go_0003712_mouse = ExTRI.update_symbols(go_0003712_mouse)

    tsv.add_field "GO:0003712" do |k,v|
      %w(human mouse rat).zip([go_0003712_human, go_0003712_mouse, go_0003712_rat])
        .select{|n,g| g.include? v.flatten.first }
        .collect{|n,g| n }
    end

    tfclass = TFClass.tfs.list
    tfclass << "AP1"
    tfclass << "NFKB"
    tsv.add_field "TFClass" do |pair,values|
      tf = pair.split(":").first
      (tfclass.include? tf) ? "TFClass" : ""
    end

    tsv.add_field "Auto-regulation" do |pair,values|
      (values[0] == values[1]) ? "Auto-regulation" : ""
    end

    tsv
  end

  dep_task :CollecTRI, ExTRI, :pairs_final

  dep :ExTRI_final
  dep :pairs_final
  task :sentence_coverage => :tsv do 
    extri = step(:ExTRI_final).path.tsv :type => :double

    extri.add_field "TF:TG" do |k,values|
      k.split(":").values_at(2,3) * ":"
    end

    extri.attach step(:pairs_final)
  end

  dep :pairs
  extension :tsv
  task :filtered_pairs => :tsv do
    tsv = step(:pairs).load
    tsv.select("TFClass" => "TFClass")

    confidence = %w(ExTRI HTRI TFactS)
    presence_fields = tsv.fields.select{|f| f.include? 'present'}

    rejected = []
    keep = []
    tsv.through do |k,v|
      present_db = presence_fields.select{|f| v[f].any? }.collect{|f| f.match(/\[(.*)\]/)[1] }

      low_conf_db = []
      confidence.each do |db|
        low_conf_db << db unless v["[#{db}] Confidence"].include? 'High'
      end

      if (present_db - low_conf_db).empty?
        rejected << k
      else
        keep << k
      end
    end
    
    tsv.select(keep)
  end

  task :tfcheckpoint_tf_class => :tsv do
    tfclass = TFClass.tfs.list
    tsv = TFCheckpoint.tfs.tsv(:merge => true)
    tsv.add_field "TFClass new" do |tf,values|
      (tfclass.include? tf) ? "TFClass new" : ""
    end
  end

  dep :ExTRI_final, :test_set => []
  dep :pairs
  task :ExTRI_coverage => :tsv do |include_HTRI|
    tsv = step(:ExTRI_final).load

    tsv.add_field "Pair" do |k,values|
      k.split(":").values_at(2,3) * ":"
    end

    pairs = step(:pairs).load
    good_fields = pairs.fields.reject{|f| f.include? "ExTRI" }

    tsv.attach pairs.slice(good_fields)
  end

end
