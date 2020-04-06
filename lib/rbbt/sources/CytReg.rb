require 'rbbt-util'
require 'rbbt/resource'

module CytReg
  extend Resource
  self.subdir = 'share/databases/CytReg'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  CytReg.claim CytReg.tf_cyt, :proc do 
    require 'rbbt/tsv/excel'
    url = "https://cytreg.bu.edu/interactions_table_full.xlsx"
    tsv = TSV.excel(url, :merge => true)
    tsv = tsv.reorder("TF",nil, :zipped => true, :merge => true)
    tsv.key_field = "Transcription Factor (Associated Gene Name)"
    tsv.fields = tsv.fields.tap{|f| f[0] = "Cytokine (Associated Gene Name)"}
    tsv.fields = tsv.fields.collect{|f| f == "PMIDs" ? "PMID" : f }
    tsv
  end
end

iif CytReg.tf_cyt.produce(true).find if __FILE__ == $0

