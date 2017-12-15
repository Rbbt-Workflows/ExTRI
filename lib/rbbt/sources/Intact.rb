require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'
require 'rbbt/sources/PRO'

module Intact
  extend Resource
  self.subdir = 'share/databases/Intact'

  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Intact.claim Intact.tf_tg, :proc do 

    uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
    gene2name = Organism.identifiers(Intact.organism).index :target => "Associated Gene Name", :order => true, :persist => true
    uni2nameHsa = UniProt.identifiers.Hsa.index :target => "Associated Gene Name", :order => true, :persist => true
    uni2nameMmu = UniProt.identifiers.Mmu.index :target => "Associated Gene Name", :order => true, :persist => true

    dumper = TSV::Dumper.new :key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "PMID", "Method ID"], :type => :double, :namespace => Intact.organism
    dumper.init
    TSV.traverse Rbbt.share.databases.FNL.Nov2017_update["Intact"], :into => dumper, :type => :line, :bar => true do |line|
      next if line =~ /^#/
      parts = line.split("\t")
      tf_ids = parts.values_at(0,2,4).collect{|str| str.split("|")}.flatten.select{|id| id.split(":").first == "uniprotkb"}.collect{|id| id.split(":").last.split("(").first}
      tg_ids = parts.values_at(1,3,5).collect{|str| str.split("|")}.flatten.select{|id| id.split(":").first == "uniprotkb"}.collect{|id| id.split(":").last.split("(").first}

      tf_ids += uni_equivalences.values_at(*tf_ids).compact.flatten
      tg_ids += uni_equivalences.values_at(*tg_ids).compact.flatten

      tf_names = tf_ids.collect{|id| gene2name[id] || uni2nameHsa[id] || uni2nameMmu[id]}.compact
      tg_names = tg_ids.collect{|id| gene2name[id] || uni2nameHsa[id] || uni2nameMmu[id]}.compact
      next if tf_names.empty? or tg_names.empty?

      pmid = parts[8].split("|").select{|f| f.include? 'pubmed'}.collect{|f| f.split(":").last} * ";"
      method = parts[6].split("|").select{|f| f.include? 'psi-mi'}.collect{|f| f.scan(/(MI:\d+)/).first} * ";"

      [tf_names.first, [tg_names.first, pmid, method]]
    end
    
    TSV.collapse_stream dumper
  end
end

if __FILE__ == $0
  require 'rbbt/workflow'
  Workflow.require_workflow "FNL"
  iif Intact.tf_tg.produce(true).find
end

