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
    tsv = tsv.reorder("TF",nil, :merge => true)
    tsv.key_field = "Transcription Factor (Associated Gene Name)"
    tsv.fields = tsv.fields.tap{|f| f[0] = "Cytokine (Associated Gene Name)"}
    tsv.fields = tsv.fields.collect{|f| f == "PMIDs" ? "PMID" : f }
    tsv
  end

  CytReg.claim CytReg.tf_cyt2, :proc do 
    require 'rbbt/tsv/csv'
    url = "https://cytreg.bu.edu/cytreg_v2.csv"
    tsv = TSV.csv(url, :merge => true).to_double
    tsv = tsv.reorder("tf",tsv.fields - %w(tf source), :one2one => true, :merge => true)
    tsv.key_field = "Transcription Factor (Associated Gene Name)"
    tsv.fields = tsv.fields.tap{|f| f[0] = "Cytokine (Associated Gene Name)"}
    tsv.fields = tsv.fields.collect{|f| f == "PMIDs" ? "PMID" : f }
    tsv
  end

  CytReg.claim CytReg.tf_cyt_merge, :proc do 
    tsv = CytReg.tf_cyt1.tsv
    tsv = tsv.merge_zip CytReg.tf_cyt2.tsv
    tsv
  end
end

Log.tsv CytReg.tf_cyt.produce(true).tsv if __FILE__ == $0
Log.tsv CytReg.tf_cyt2.produce(true).tsv if __FILE__ == $0
Log.tsv CytReg.tf_cyt_merge.produce(true).tsv if __FILE__ == $0

