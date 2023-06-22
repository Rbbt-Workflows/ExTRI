require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'
require 'rbbt/sources/PRO'

module IntAct
  extend Resource
  self.subdir = 'share/databases/IntAct'

  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  IntAct.claim IntAct.tf_tg, :proc do 

    uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
    gene2name = Organism.identifiers(IntAct.organism).index :target => "Associated Gene Name", :order => true, :persist => true

    gene2uniHsa = Organism.identifiers(IntAct.organism("Hsa")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniMmu = Organism.identifiers(IntAct.organism("Mmu")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniRno = Organism.identifiers(IntAct.organism("Rno")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true

    uni2nameHsa = UniProt.identifiers.Hsa.index :target => "Associated Gene Name", :order => true, :persist => true
    #uni2nameMmu = UniProt.identifiers.Mmu.index :target => "Associated Gene Name", :order => true, :persist => true
    #uni2nameRno = UniProt.identifiers.Rno.index :target => "Associated Gene Name", :order => true, :persist => true

    dumper = TSV::Dumper.new :key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "PMID", "Method ID"], :type => :double, :namespace => IntAct.organism
    dumper.init
    parser = TSV::Parser.new(Rbbt.share.databases.ExTRI.Nov2017_update["Intact"].find)
    fields = parser.fields
    TSV.traverse parser, :into => dumper, :type => :line, :bar => true do |line|
      next if line =~ /^#/
      parts = line.split("\t")

      tf_ids = parts.values_at(0,2,4).collect{|str| str.split("|")}.flatten
      tf_ids = tf_ids.select{|id| id.split(":").first == "uniprotkb"}.collect{|id| id.split(":").last.split("(").first}.uniq
      tf_ids += tf_ids.collect{|id| gene2uniHsa[id] || gene2uniMmu[id] || gene2uniRno[id] }.compact.uniq
      tf_ids += uni_equivalences.values_at(*tf_ids).compact.flatten.uniq
      tf_names = gene2name.values_at(*tf_ids).compact.flatten.uniq

      tg_ids = parts.values_at(1,3,5).collect{|str| str.split("|")}.flatten
      tg_ids = tg_ids.select{|id| id.split(":").first == "uniprotkb"}.collect{|id| id.split(":").last.split("(").first}.uniq
      tg_ids += tg_ids.collect{|id| gene2uniHsa[id] || gene2uniMmu[id] || gene2uniRno[id] }.compact.uniq
      tg_ids += uni_equivalences.values_at(*tg_ids).compact.flatten.uniq
      tg_names = gene2name.values_at(*tg_ids).compact.flatten.uniq

      pmid = parts[8].split("|").select{|f| f.include? 'pubmed'}.collect{|f| f.split(":").last} * ";"

      next if tf_names.empty? or tg_names.empty?

      type_a, type_b = parts.values_at 20, 21

      if type_a.include?('gene') && type_b.include?('protein')
        tf_names, tg_names = tg_names, tf_names
      end

      method = parts[6].split("|").select{|f| f.include? 'psi-mi'}.collect{|f| f.scan(/(MI:\d+)/).first} * ";"

      [tf_names.first, [tg_names.first, pmid, method]]
    end
    
    TSV.collapse_stream dumper
  end
end

if __FILE__ == $0
  require 'rbbt/workflow'
  Workflow.require_workflow "ExTRI"
  iif IntAct.tf_tg.produce(true).find
end

