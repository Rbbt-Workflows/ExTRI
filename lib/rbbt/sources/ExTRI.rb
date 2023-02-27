require 'rbbt-util'
require 'rbbt/resource'

require 'rbbt/sources/organism'
require 'rbbt/sources/PRO'

module ExTRI
  extend Resource
  self.subdir = 'share/databases/ExTRI'
  self.set_libdir
  self.search_paths.merge!(:default => :lib)

  def self.organism
    "Hsa/oct2016"
  end

  def self.updated_organism
    "Hsa/feb2023"
  end


  def self.entrez_gene_index

    entrez_gene_index = {}
    TSV.traverse CMD.cmd('grep "^\(9606\|10090\|10116\)"', :in => Rbbt.share.databases.entrez.gene_info.open, :pipe => true), :type => :array do |line|
      parts = line.split("\t", -1)
      entrez = parts[1]
      external = parts[5]
      codes = external.split("|")

      codes.each do |code|
        code = code.sub(/Ensembl:/,'')
        entrez_gene_index[code] = entrez
      end
    end

    entrez_gene_index
  end


  ExTRI.claim ExTRI.GOA, :proc do 
    #url = "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_ref_uniprot.gz"
    #text = CMD.cmd('grep has_regulation_target | grep "taxon:\(9606\|10090\|10116\)" | grep $\'\t\(EXP\|IDA\|IPI\|IMP\|IGI\|IEP\|TAS\)\t\'', :in => Open.open(url), :pipe => true)
    #taxons = [9606, 10090, 10116]
    #entrez_gene_index = Rbbt.share.databases.entrez.gene_info.index :target => [1], :fields => [5], :grep => /^(#{taxons * "|"})/

    entrez_gene_index = ExTRI.entrez_gene_index

    TSV.traverse CMD.cmd('grep "^\(9606\|10090\|10116\)"', :in => Rbbt.share.databases.entrez.gene_info.open, :pipe => true), :type => :array do |line|
      parts = line.split("\t", -1)
      entrez = parts[1]
      external = parts[5]
      codes = external.split("|")

      codes.each do |code|
        entrez_gene_index[code] = entrez
      end
    end

    text = Rbbt.data['.source']['gene_association.goa_ref_uniprot.filtered'].open

    mouse_index = Organism.identifiers("Mmu").index(:target => "Entrez Gene ID", :persist => true, :order => true)
    rat_index = Organism.identifiers("Rno").index(:target => "Entrez Gene ID", :persist => true, :order => true)
    human_index = Organism.identifiers(ExTRI.organism("Hsa")).index(:target => "Entrez Gene ID", :persist => true, :order => true)


    tsv = TSV.setup({}, :key_field => "UniProt/SwissProt Accession", :fields => ["Entrez Gene ID", "GO ID", "Evidence type", "Provider"], :type => :double)

    TSV.traverse text, :type => :array, :into => tsv do |line|
      parts = line.split("\t", -1)
      uniprot = parts[1]
      go_term = parts[4]
      evidence = parts[6]
      provider = parts[14]
      interaction_info = parts[15]
      raw_targets = interaction_info.scan(/has_regulation_target\((.*?)\)/).flatten
      targets = raw_targets.collect{|rt|
        ns,_sep, code = rt.partition ":"
        case ns
        when "MGI"
          entrez_gene_index[code] || mouse_index[code] || "MISSING(#{code})"
        when "UniProtKB"
          "UNIPROT(#{code})"
        else
          entrez_gene_index[code] ||   human_index[code] || mouse_index[code] || rat_index[code] || "MISSING(#{code})"
        end
      }.compact.reject{|t| t == ""}

      next if targets.empty?

      num_targets = targets.length

      [uniprot, [targets, [go_term]*num_targets, [evidence]*num_targets, [provider]*num_targets]]
    end

    tsv
  end

  ExTRI.claim ExTRI.HTRIDB, :proc do 
    file = Rbbt.data['.source']['raw_data_14-08-04.csv']
    tsv = TSV.open(file, :header_hash => '', :sep => ';', :fix => Proc.new{|l| l[0..-1]}, :key_field => "GENEID_TF", :fields => "GENEID_TG;TECHNIQUE;PUBMED_ID".split(";"), :merge => true)
    tsv.key_field = "Entrez Gene ID"
    tsv.fields = ["Entrez Gene ID", "Technique", "PMID"]


    tsv.namespace = Organism.default_code('Hsa')
    tsv.identifiers = Organism.identifiers(tsv.namespace)

    tsv = tsv.change_key("UniProt/SwissProt Accession")

    tsv
  end

  ExTRI.claim ExTRI.HTRIDB_HT_methods, :proc do 
    all = ExTRI.HTRIDB.tsv.column("Technique").values.flatten.uniq
    small = Rbbt.data['.source']['small_scale_experiments_only.csv'].tsv(:header_hash => '', :sep => ';', :merge => true).column("TECHNIQUE").values.flatten.compact.uniq

    (all - small) * "\n"
  end

  ExTRI.claim ExTRI.Encode, :proc do 
    url = 'http://encodenets.gersteinlab.org/enets2.Proximal_filtered.txt'
    tsv = TSV.open(url, :fields => [2,1], :merge => true, :sep => /\s+/)
    tsv.key_field = "Associated Gene Name"
    tsv.fields = ["Associated Gene Name", "Type"]

    url2 = 'http://encodenets.gersteinlab.org/enets3.Distal.txt'
    tsv2 = TSV.open(url2, :fields => [2,1], :merge => true, :sep => /\s+/)
    tsv2.key_field = "Associated Gene Name"
    tsv2.fields = ["Associated Gene Name", "Type"]

    tsv2.each do |k,v|
      if old = tsv[k]
        old = old.dup
        old[0].concat v[0]
        old[1].concat v[1]
        tsv[k] = old
      else
        tsv[k] = v
      end
    end

    tsv.namespace = Organism.default_code('Hsa')
    tsv.identifiers = Organism.identifiers(tsv.namespace)
    tsv = tsv.change_key("UniProt/SwissProt Accession")
    tsv = tsv.swap_id("Associated Gene Name", "Entrez Gene ID")

    tsv
  end

  ExTRI.claim ExTRI.IntAct, :proc do
    mouse_index = Organism.identifiers("Mmu").index(:target => "Entrez Gene ID", :persist => true, :order => true)
    rat_index = Organism.identifiers("Rno").index(:target => "Entrez Gene ID", :persist => true, :order => true)
    human_index = Organism.identifiers(ExTRI.organism("Hsa")).index(:target => "Entrez Gene ID", :persist => true, :order => true)

    entrez_gene_index = ExTRI.entrez_gene_index

    tsv = TSV.setup({}, :key_field => "UniProt/SwissProt Accession", :fields => ["Entrez Gene ID", "Method ID", "PMID"], :type => :double)
    TSV.traverse Rbbt.data['.source']['mammalian_protein_gene_interactions.intact_14-09-03.txt'], :type => :array, :into => tsv do |line|
      uni,ensembl,method,pmid = line.split("\t")
      uni = uni.split(":").last
      ensembl = ensembl.split(":").last
      entrez = entrez_gene_index[ensembl] || human_index[ensembl] || mouse_index[ensembl] || rat_index[ensembl] || "MISSING(#{ensembl}"
      pmid = pmid.split(":").last
      method = method.scan(/"(.*)"/).flatten.first
      [uni, [[entrez], [method], [pmid]]]
    end

    tsv
  end

  ExTRI.claim ExTRI.TFACTS, :proc do
    tsv = TSV.setup({}, :key_field => "UniProt/SwissProt Accession", :fields => ["Entrez Gene ID", "Provider", "PMID"], :type => :double)

    mouse_index_uni = Organism.identifiers("Mmu").index(:target => "UniProt/SwissProt Accession", :persist => true, :order => true)
    rat_index_uni = Organism.identifiers("Rno").index(:target => "UniProt/SwissProt Accession", :persist => true, :order => true)
    human_index_uni = Organism.identifiers(ExTRI.organism("Hsa")).index(:target => "UniProt/SwissProt Accession", :persist => true, :order => true)

    TSV.traverse Rbbt.data['.source']['TFactS_2012-08-09.trancritption_factor_entrez-target_gene_entrez-provenance-organism-pmid_ref.txt'], :type => :array, :into => tsv do |line|
      next if line =~ /^#/
      tg,tf,provider,organism,pmid = line.split("\t")

      pmid = "" if pmid.nil?
      provider = provider.gsub(/\s/,'').gsub(/^;/,'').gsub(/;$/,'')
      pmid = pmid.gsub(/\s/,'').gsub(/^;/,'').gsub(/;$/,'')
      tf = human_index_uni[tf] || mouse_index_uni[tf] || rat_index_uni[tf] || "MISSING(#{tf})"
      [tf,[ [tg], [provider], [pmid]]]
    end

    tsv
  end

  ExTRI.claim ExTRI.TFactS_flagged_articles, :proc do
    threshold = 10
    pmid_counts = {}
    TSV.traverse TFactS.tf_tg do |tf,values|
      Misc.zip_fields(values).each do |target,sign,species,source,pmids|
        pair = [tf,target] * ":"
        next if pmids.nil? or pmids.empty?
        pmids.split(";").each do |pmid|
          pmid_counts[pmid] ||= []
          pmid_counts[pmid] << pair
        end
      end
    end

    flagged = pmid_counts.select{|pmid,list| list.uniq.length > threshold}.collect{|pmid,list| pmid}

    flagged * "\n"
  end

  ExTRI.claim ExTRI.Thomas2015, :proc do
    uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
    mgi2uni = Organism.identifiers("Mmu/feb2014").index :target => "Associated Gene Name", :persist => true
    uni2name = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :persist => true

    CaseInsensitiveHash.setup(mgi2uni)
    CaseInsensitiveHash.setup(uni2name)

    thomas = ExTRI.Nov2017_update.Thomas2015.tsv(:key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "sentence", "class", "details", "PMID"], :merge => true, :type => :double)

    new = thomas.annotate({})

    thomas.through do |tf,values|
      Misc.zip_fields(values).each do |tg,*rest|

        tf_alts = [tf, mgi2uni[tf]].compact
        tf_alts += uni_equivalences.values_at(*tf_alts).flatten.compact
        tf_new = uni2name.values_at(*tf_alts).compact.first

        tg_alts = [tg, mgi2uni[tg]].compact.flatten
        tg_alts += uni_equivalences.values_at(*tg_alts).flatten.compact
        tg_new = uni2name.values_at(*tg_alts).compact.first

        next if tg_new.nil? or tf_new.nil?
        new_values = [tg_new] + rest
        new.zip_new(tf_new, new_values)
      end
    end

    new
  end

  ExTRI.claim ExTRI.NTNU_curated, :proc do
    dumper = TSV::Dumper.new :key_field => "TRI", :fields => ["Transcription Factor (Associated Gene Name)", "Target Gene (Associated Gene Name)", "Sign", "PMID"], :type => :list
    dumper.init
    TSV.traverse Rbbt.data["ExTRI_curated.tsv"], :type => :list, :into => dumper do |tri,values|
      valid, sign, negated, *rest = values
      pmid, number, tf, tg = tri.split(":")

      sign = case sign
             when "+"
               "UP"
             when "-"
               "DOWN"
             else
               ""
             end

      next if negated == "true" && sign == ""

      sign = "" if negated == "true"

      [tri, [tf, tg, sign, pmid]]
    end
  end

end


