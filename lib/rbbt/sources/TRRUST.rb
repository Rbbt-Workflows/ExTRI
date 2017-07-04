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


  TRRUST.claim TRRUST.tf_tg, :proc do 
    url = "http://www.grnpedia.org/trrust/trrust_rawdata.txt"
    tsv = TSV.open(url, :merge => true)
    tsv.key_field = "Transcription Factor (Associated Gene Name)"
    tsv.fields = ["Target Gene (Associated Gene Name)", "Regulation", "PMID"]
    tsv.namespace = TRRUST.organism
    tsv.to_s
  end
end

iif TRRUST.tf_tg.produce.find if __FILE__ == $0

