require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'

module TFCheckpoint
  extend Resource
  self.subdir = 'share/databases/TFCheckpoint'

  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  TFCheckpoint.claim TFCheckpoint.tfs, :proc do 
    url = "http://www.tfcheckpoint.org/data/TFCheckpoint_download_180515.txt"
    tsv = TSV.open(url, :sep2 => /[,|]\s*/, :header_hash => "")
    tsv.key_field = "Associated Gene Name"
    tsv.namespace = TFCheckpoint.organism
    tsv.to_s
  end
end

iif TFCheckpoint.tfs.produce(true).find if __FILE__ == $0

