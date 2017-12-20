module FNL

  PSI_MITAB_FIELDS =<<-EOF.split("\n").collect{|l| l.strip}
    Unique identifier for interactor A
    Unique identifier for interactor B
    Alternative identifier for interactor A
    Alternative identifier for interactor B
    Aliases for A
    Aliases for B
    Interaction detection methods
    First author
    Identifier of the publication
    NCBI Taxonomy identifier for interactor A
    NCBI Taxonomy identifier for interactor B
    Interaction types
    Source databases and identifiers
    Interaction identifier(s)
    Confidence score
    Complex expansion
    Biological role A 
    Biological role B
    Experimental role A 
    Experimental role B 
    Interactor type A 
    Interactor type B
    Xref for interactor A
    Xref for interactor B
    Xref for the interaction
    Annotations for interactor A
    Annotations for Interactor B
    Annotations for the interaction
    NCBI Taxonomy identifier for the host organism
    Parameters of the interaction
    Creation date
    Update date
    Checksum for interactor A
    Checksum for interactor B
    Checksum for interaction
    negative
    Feature(s) for interactor A
    Feature(s) for interactor B
    Stoichiometry for interactor A
    Stoichiometry for interactor B
    Participant identification method for interactor A
    Participant identification method for interactor B
  EOF

  
  dep :sentence_coverage_NER
  task :PSI_MITAB_FNL => :tsv do
    psi_key_field, *psi_fields = PSI_MITAB_FIELDS

    dumper = TSV::Dumper.new(:key_field => psi_key_field, :fields => psi_fields, :type => :double)
    dumper.init

    name2ens = Organism.identifiers(FNL.organism).index :target => "Ensembl Gene ID", :order => true, :persist => true
    name2uni = Organism.identifiers(FNL.organism).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    name2ent = Organism.identifiers(FNL.organism).index :target => "Entrez Gene ID", :order => true, :persist => true
    ens2enst = Organism.transcripts(FNL.organism).tsv :key_field => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true, :type => :single
    enst2ensp = Organism.transcripts(FNL.organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Protein ID"], :persist => true, :type => :single
    ensg2penst = Appris.principal_isoforms.tsv :fields => ["Ensembl Transcript ID"], :type => :flat, :persist => true
    ensp2uni = Organism.protein_identifiers(FNL.organism).index :target => "UniProt/SwissProt Accession", :persist => true

    parser = TSV::Parser.new step(:sentence_coverage_NER)
    sent_key_field, *sent_fields = parser.all_fields
    TSV.traverse parser, :into => dumper do |sent_id, sent_values|
      values = {}

      id = "ntnu:" << sent_id.gsub(":", "_")

      confidence = sent_values[sent_fields.index("Prediction confidence")]
      sentence = sent_values[sent_fields.index("Sentence")]

      pmid, sent, tf_name, tg_name = sent_id.split(":")

      NamedArray.setup(sent_values, sent_fields)

      tf_tax = sent_values["TF Tax"] 
      tg_tax = sent_values["TG Tax"]
      tf_tax = "9606" if tf_tax.nil? or tf_tax.empty?
      tg_tax = "9606" if tg_tax.nil? or tg_tax.empty?

      pref_tax = "taxid:"
      pref_uni = "uniprotkb:"
      pref_ens = "ensembl:"
      pref_ent = "entrezgene/locuslink:"
      pref_name = "ensembl:"

      tf_ens, tg_ens = name2ens.values_at tf_name, tg_name
      tf_ent, tg_ent = name2ent.values_at tf_name, tg_name
      tf_uni, tg_uni = name2uni.values_at tf_name, tg_name

      next if tf_ens.nil?

      tf_enst = (ensg2penst[tf_ens] || []).first || ens2enst[tf_ens]
      tf_enp = enst2ensp[tf_enst] 
      tf_enp ||= enst2ensp[ens2enst[tf_ens]]

      tf_uni = ensp2uni[tf_enp] if ensp2uni[tf_enp] and ensg2penst[tf_ens]

      tf_uni ||= "Missing Uniprot ID for #{tf_name}"
      tg_uni ||= "Missing Uniprot ID for #{tg_name}"

      tf_ens ||= "Missing Ensembl Gene ID for #{tf_name}"
      tg_ens ||= "Missing Ensembl Gene ID for #{tg_name}"

      tf_ent ||= "Missing Entrez Gene ID for #{tf_name}"
      tg_ent ||= "Missing Entrez Gene ID for #{tg_name}"

      values["Unique identifier for interactor B"] = [pref_ens + tg_ens]
      values["Alternative identifier for interactor A"] = [pref_ens + tf_enp, pref_ent + tf_ent]
      values["Alternative identifier for interactor B"] = [pref_ent + tf_ent]
      values["Aliases for A"] = [pref_name + tf_name + "(gene name)"] 
      values["Aliases for B"] = [pref_name + tg_name + "(gene name)"]
      values["Interaction detection methods"] = ['psi-mi:"MI:0087"(predictive tm)']
      values["First author"] = nil
      values["Identifier of the publication"] = ["pmid:" << pmid]
      values["NCBI Taxonomy identifier for interactor A"] = ["taxid:" << tf_tax]
      values["NCBI Taxonomy identifier for interactor B"] = ["taxid:" << tg_tax]
      values["Interaction types"] = nil
      values["Source databases and identifiers"] = ['psi-mi:"MI:1264"(ntnu)']
      values["Interaction identifier(s)"] = [id]
      values["Confidence score"] = ["confidence:" + confidence]
      values["Complex expansion"] = nil
      values["Biological role A "] = ['psi-mi:"MI:2274"(regulator)']
      values["Biological role B"] =  ['psi-mi:"MI:2275"(regulator target)']
      values["Experimental role A "] = ['psi-mi:"MI:0499"(unpecified role)'] 
      values["Experimental role B "] = ['psi-mi:"MI:0499"(unpecified role)']
      values["Interactor type A "] = ['psi-mi:"MI:0326"(protein)']
      values["Interactor type B"] =  ['psi-mi:"MI:0250"(gene)']
      values["Xref for interactor A"] = nil
      values["Xref for interactor B"] = nil
      values["Xref for the interaction"] = nil
      values["Annotations for interactor A"] = nil
      values["Annotations for Interactor B"] = nil
      values["Annotations for the interaction"] = ['comment:"' + sentence + '"']
      values["NCBI Taxonomy identifier for the host organism"] = ["taxid:" << tf_tax]
      values["Parameters of the interaction"] = nil
      values["Creation date"] = nil
      values["Update date"] = nil
      values["Checksum for interactor A"] = nil
      values["Checksum for interactor B"] = nil
      values["Checksum for interaction"] = nil
      values["negative"] = nil
      values["Feature(s) for interactor A"] = nil
      values["Feature(s) for interactor B"] = nil
      values["Stoichiometry for interactor A"] = nil
      values["Stoichiometry for interactor B"] = nil
      values["Participant identification method for interactor A"] = nil
      values["Participant identification method for interactor B"] = nil

      [pref_uni + tf_uni, 
       psi_fields.collect{|f| (values[f] || ["-"]) * "|" } 
      ]
    end
  end

  #{{{ ALL PAIRS

  
  #dep :pairs
  #task :PSI_MITAB => :tsv do
  #  psi_key_field, *psi_fields = PSI_MITAB_FIELDS

  #  dumper = TSV::Dumper.new(:key_field => psi_key_field, :fields => psi_fields, :type => :double)
  #  dumper.init

  #  name2ens = Organism.identifiers(FNL.organism).index :target => "Ensembl Gene ID", :order => true, :persist => true
  #  name2uni = Organism.identifiers(FNL.organism).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
  #  name2ent = Organism.identifiers(FNL.organism).index :target => "Entrez Gene ID", :order => true, :persist => true

  #  parser = TSV::Parser.new step(:pairs)
  #  pairs_key_field, *pairs_fields = parser.all_fields
  #  TSV.traverse parser, :into => dumper do |pairs_id, pairs_values|
  #    values = {}

  #    id = "FNL:" << pairs_id.gsub(":", "_")

  #    tf_name, tg_name = pairs_id.split(":")

  #    NamedArray.setup(pairs_values, pairs_fields)

  #    pmids = pairs_values.values_at(*pairs_pmid_fields)

  #    tf_tax = pairs_values["TF Tax"] 
  #    tg_tax = pairs_values["TG Tax"]
  #    tf_tax = "9606" if tf_tax.nil? or tf_tax.empty?
  #    tg_tax = "9606" if tg_tax.nil? or tg_tax.empty?

  #    pref_tax = "taxid:"
  #    pref_uni = "uniprot/swiss-prot:"
  #    pref_ens = "ensembl:"
  #    pref_ent = "entrezgene/locuslink:"

  #    tf_ens, tg_ens = name2ens.values_at tf_name, tg_name
  #    tf_uni, tg_uni = name2uni.values_at tf_name, tg_name
  #    tf_ent, tg_ent = name2ent.values_at tf_name, tg_name

  #    tf_uni ||= "Missing Uniprot ID for #{tf_name}"
  #    tg_uni ||= "Missing Uniprot ID for #{tg_name}"

  #    tf_ens ||= "Missing Ensembl Gene ID for #{tf_name}"
  #    tg_ens ||= "Missing Ensembl Gene ID for #{tg_name}"

  #    tf_ent ||= "Missing Entrez Gene ID for #{tf_name}"
  #    tg_ent ||= "Missing Entrez Gene ID for #{tg_name}"



  #    values["Unique identifier for interactor A"] = [pref_uni + tf_uni]
  #    values["Unique identifier for interactor B"] = [pref_ent + tf_ent]
  #    values["Alternative identifier for interactor A"] = [pref_uni + tf_uni]
  #    values["Alternative identifier for interactor B"] = nil
  #    values["Aliases for A"] = nil
  #    values["Aliases for B"] = nil
  #    values["Interaction detection methods"] = nil
  #    values["First author"] = nil
  #    values["Identifier of the publication"] = ["pmid:" << pmid]
  #    values["NCBI Taxonomy identifier for interactor A"] = ["taxid:" << tf_tax]
  #    values["NCBI Taxonomy identifier for interactor B"] = ["taxid:" << tg_tax]
  #    values["Interaction types"] = nil
  #    values["Source databases and identifiers"] = nil
  #    values["Interaction identifier(s)"] = nil
  #    values["Confidence score"] = nil
  #    values["Complex expansion"] = nil
  #    values["Biological role A "] = nil
  #    values["Biological role B"] = nil
  #    values["Experimental role A "] = nil
  #    values["Experimental role B "] = nil
  #    values["Interactor type A "] = nil
  #    values["Interactor type B"] = nil
  #    values["Xref for interactor A"] = nil
  #    values["Xref for interactor B"] = nil
  #    values["Xref for the interaction"] = nil
  #    values["Annotations for interactor A"] = nil
  #    values["Annotations for Interactor B"] = nil
  #    values["Annotations for the interaction"] = nil
  #    values["NCBI Taxonomy identifier for the host organism"] = ["taxid:" << tf_tax]
  #    values["Parameters of the interaction"] = nil
  #    values["Creation date"] = nil
  #    values["Update date"] = nil
  #    values["Checksum for interactor A"] = nil
  #    values["Checksum for interactor B"] = nil
  #    values["Checksum for interaction"] = nil
  #    values["negative"] = nil
  #    values["Feature(s) for interactor A"] = nil
  #    values["Feature(s) for interactor B"] = nil
  #    values["Stoichiometry for interactor A"] = nil
  #    values["Stoichiometry for interactor B"] = nil
  #    values["Participant identification method for interactor A"] = nil
  #    values["Participant identification method for interactor B"] = nil

  #    [tf_ens, psi_fields.collect{|f| (values[f] || ["-"]) * "|" } ]
  #  end
  #end
end

