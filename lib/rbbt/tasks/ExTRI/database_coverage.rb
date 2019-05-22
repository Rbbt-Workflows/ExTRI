require 'rbbt/sources/TFClass'
module ExTRI
  AP1_SYN=%w(FOS FOSB JUN JUNB JUND FOSL1 FOSL2)
  NFKB_SYN=%w(NFKB1 NFKB2 RELA RELB)


  DATABASES=%w(ExTRI HTRI TFacts TRRUST Intact Encode)

  helper :normalize_db do |db|
    new = db.annotate({})
    translation = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :order => true, :persist => true
    
    TSV.traverse db, :into =>  new do |pair,values|
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
        extra = [name] + Misc.zip_fields(all_values).collect{|l| l*";"} 
      else
        extra = [""] + [""] * num_new_fields
      end

      new_values = values + extra

      [k, new_values]
    end
  end

  dep :ExTRI_confidence, :test_set => []
  task :sentence_coverage => :tsv do
    id_file = Organism.identifiers(ExTRI.organism)

    encode = ExTRI.Encode.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip

    #goa = ExTRI.GOA.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
    goa = GO.tf_tg.tsv(:merge => true).unzip

    #intact = ExTRI.Intact.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
    intact = Intact.tf_tg.tsv(:merge => true).unzip

    tfacts = TFacts.tf_tg.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :merge => true, :zipped => true).unzip
    trrust = TRRUST.Hsa.tf_tg.tsv(:merge => true).unzip
    htri = HTRI.tf_tg.tsv(:merge => true).unzip
    signor = Signor.tf_tg.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => UniProt.identifiers.Hsa).unzip
    thomas = ExTRI.Thomas2015.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "sentence", "class", "details", "PMID"], :merge => true).unzip

    flagged = ExTRI.TFacts_flagged_articles.list
    tfacts.add_field "Confidence" do |tf,values|
      sign,species,source,pmids = values
      (source.downcase == "pubmed" and (pmids.split(";") - flagged).empty?) ? "Low" : "High"
    end
  
    tsv = step(:ExTRI_confidence).load

    tsv = attach_db tsv, htri, "HTRI"
    tsv = attach_db tsv, trrust, "TRRUST"
    tsv = attach_db tsv, tfacts, "TFacts"
    #tsv = attach_db tsv, encode, "Encode"
    tsv = attach_db tsv, goa, "GOA"
    tsv = attach_db tsv, intact, "Intact"
    tsv = attach_db tsv, signor, "Signor"
    #tsv = attach_db tsv, thomas, "Thomas2015"

    tsv
  end

  desc <<-EOF 

List all TF:TG pairs across ExTRI and other resources along with confidence estimates and other information from those resources.

The confidence estimate for ExTRI pairs uses by default 2 PMIDs or 2 sentences or a score over 1.6.

  EOF
  dep :ExTRI_confidence, :pmids => 2, :sentences => 2, :score => 1.6, :test_set => []
  task :pairs => :tsv do
    id_file = Organism.identifiers(ExTRI.organism)

    encode = ExTRI.Encode.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip

    #goa = ExTRI.GOA.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
    goa = GO.tf_tg.tsv(:merge => true).unzip

    #intact = ExTRI.Intact.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => id_file).swap_id("Entrez Gene ID", "Associated Gene Name", :identifiers => id_file).unzip
    intact = Intact.tf_tg.tsv(:merge => true).unzip

    tfacts = TFacts.tf_tg.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :merge => true, :zipped => true).unzip
    trrust = TRRUST.Hsa.tf_tg.tsv(:merge => true).unzip
    htri = HTRI.tf_tg.tsv(:merge => true).unzip
    signor = Signor.tf_tg.tsv(:merge => true).change_key("Associated Gene Name", :identifiers => UniProt.identifiers.Hsa).unzip
    #thomas = ExTRI.Thomas2015.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "class", "details", "sentence", "PMID"], :merge => true).unzip
    #cp = TFCheckpoint.tfs.tsv(:merge => true)

    flagged = ExTRI.TFacts_flagged_articles.list
    tfacts.add_field "Confidence" do |tf,values|
      sign,species,source,pmids = values
      (source.downcase == "pubmed" and (pmids.split(";") - flagged).empty?) ? "Low" : "High"
    end

    orig = step(:ExTRI_confidence).load

    tsv = TSV.setup({}, :key_field => "TF:TG", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "[ExTRI] Confidence", "[ExTRI] PMID"], :type => :list, :namespace => ExTRI.organism)

    confidence = orig.fields.select{|f| f.include? "Thresh"}.first
    pmids = {}
    conf = {}
    orig.through do |key,values|
      c = values[confidence]
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
      [tfacts, "TFacts"],
      #[encode, "Encode"], # Don't use Encode
      [goa, "GOA"],
      [intact, "Intact"],
      [signor, "Signor"],
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

    tsv.process "Transcription Factor (Associated Gene Name)" do |list,key,values|
      key.split(":").first
    end

    tsv.process "Target Gene (Associated Gene Name)" do |list,key,values|
      key.split(":").last
    end

    #cp.add_field "present" do
    #  "TFCheckpoint"
    #end

    #cp.fields = cp.fields.collect{|f| "[TFCheckpoint] " << f}
    #cp.key_field = "TF"
    #
    #new_cp = cp.annotate({})

    #translation = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :order => true, :persist => true
    #TSV.traverse cp, :into => new_cp do |tf, values|
    #  new_tf = translation[tf] 
    #  if new_tf.nil?
    #    Log.error "Gene name not found from TFCheckpoint: " + tf
    #    new_tf = tf
    #  end
    #  Log.info "TFCheckpoint Translation " << [tf, new_tf] * " => " if tf != new_tf
    #  [new_tf, values]
    #end

    #cp = new_cp

    #tsv.attach cp

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

  task :tfcheckpoint_tf_class => :tsv do
    tfclass = TFClass.tfs.list
    tsv = TFCheckpoint.tfs.tsv(:merge => true)
    tsv.add_field "TFClass new" do |tf,values|
      (tfclass.include? tf) ? "TFClass new" : ""
    end
  end


end
