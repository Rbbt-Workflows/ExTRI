require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'

module TRRUST
  extend Resource
  self.subdir = 'share/databases/TRRUST'

  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  TRRUST.claim TRRUST.Hsa.tf_tg, :proc do 
    url = "http://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
    tsv = TSV.open(url, :merge => true)
    tsv.key_field = "Transcription Factor (Associated Gene Name)"
    tsv.fields = ["Target Gene (Associated Gene Name)", "Regulation", "PMID"]
    tsv.namespace = TRRUST.organism
    tsv.to_s
  end

  TRRUST.claim TRRUST.Mmu.tf_tg, :proc do 
    url = "http://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv"
    tsv = TSV.open(url, :merge => true)
    tsv.key_field = "Transcription Factor (Associated Gene Name)"
    tsv.fields = ["Target Gene (Associated Gene Name)", "Regulation", "PMID"]
    tsv.namespace = TRRUST.organism
    tsv.to_s
  end
end

iif TRRUST.Hsa.tf_tg.produce.find if __FILE__ == $0
iif TRRUST.Mmu.tf_tg.produce.find if __FILE__ == $0

