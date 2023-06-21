require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'

module Pavlidis
  extend Resource
  self.subdir = 'share/databases/Pavlidis'

  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Pavlidis.claim Pavlidis.tf_tg_OLD, :proc do 
    url = "https://doi.org/10.1371/journal.pcbi.1009484.s034"

    ent2ensg_hsa = Organism.identifiers(organism("Hsa")).index :target => "Ensembl Gene ID", :fields => ["Entrez Gene ID"], :persist => true
    ent2ensg_mmu = Organism.identifiers(organism("Mmu")).index :target => "Ensembl Gene ID", :fields => ["Entrez Gene ID"], :persist => true
    ent2ensg_rno = Organism.identifiers(organism("Rno")).index :target => "Ensembl Gene ID", :fields => ["Entrez Gene ID"], :persist => true

    ensg2name = Organism.identifiers(organism("Hsa")).index :target => "Associated Gene Name", :fields => ["Ensembl Gene ID"], :persist => true

    mmu2hsa = Organism.ortholog_Hsa(organism("Mmu")).tsv(:type => :flat)
    rno2hsa = Organism.ortholog_Hsa(organism("Mmu")).tsv(:type => :flat)

    key_field, *fields =<<-EOF.split("\n")
Transcription Factor (Associated Gene Name)
Target Gene (Associated Gene Name)
Experiment Type
Experimental Method
Mode of action
Curation
    EOF
    tsv = TSV.setup({}, :key_field => key_field, :fields => fields, :type => :double)

    TSV.traverse Open.open(url), :type => :array, :into => tsv do |line|
      next unless line =~ /^tri_/

      id, exp_id, tf, tg, exp_type, int, int_mi, tf_role, tf_role_mi,
        target_role, target_role_mi, exp_method, exp_method_mi, mode, mode_mi,
        curation, curation_mi = line.split("\t")

      tf_ens = ent2ensg_hsa[tf] || ent2ensg_mmu[tf] || ent2ensg_rno[tf]
      tg_ens = ent2ensg_hsa[tg] || ent2ensg_mmu[tg] || ent2ensg_rno[tg]

      if tf_ens.nil? || tg_ens.nil?
        Log.warn ["Ens", tf, tf_ens, tg, tg_ens] * ", "
        next
      end

      tf_hsa = mmu2hsa[tf_ens] || rno2hsa[tf_ens] || tf_ens
      tg_hsa = mmu2hsa[tg_ens] || rno2hsa[tg_ens] || tg_ens

      tf_hsa = tf_hsa.first if Array === tf_hsa
      tg_hsa = tg_hsa.first if Array === tg_hsa

      if tf_hsa.nil? || tg_hsa.nil?
        Log.warn ["Hsa", tf, tf_ens, tf_hsa, tg, tg_ens, tg_hsa] * ", "
        next
      end

      tf_name = ensg2name[tf_hsa]
      tg_name = ensg2name[tg_hsa]

      if tf_name.nil? || tg_name.nil?
        Log.warn ["Name", tf, tf_ens, tf_hsa, tf_name, tg, tg_ens, tg_hsa, tg_name] * ", "
        next
      end

      values = [[tg_name], [exp_type], [exp_method], [mode]]

      [tf_name, values] 
    end

    tsv
  end

  Pavlidis.claim Pavlidis.tf_tg, :proc do 
    url = "https://doi.org/10.1371/journal.pcbi.1009484.s025"

    key_field, *fields =<<-EOF.split("\n")
Transcription Factor (Associated Gene Name)
Target Gene (Associated Gene Name)
PMID
Mode of action
    EOF
    tsv = TSV.setup({}, :key_field => key_field, :fields => fields, :type => :double)

    TSV.traverse Open.open(url), :type => :array, :into => tsv do |line|
      next unless line =~ /^tri_/
      parts = line.split("\t")

      tf = parts[2]
      tg = parts[3]
      pmid = parts[8]
      mode = parts[15]

      values = [tg, pmid, mode]

      [tf, values]
    end

    tsv
  end

end
