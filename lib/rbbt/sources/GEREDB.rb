require 'rbbt-util'
require 'rbbt/resource'
require_relative 'TFClass'

module GEREDB
  extend Resource
  self.subdir = 'share/databases/GEREDB'

  def self.organism(org="Hsa")
    require 'rbbt/sources/organism'
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  GEREDB.claim GEREDB[".source"].data, :proc do 
    tsv = TmpFile.with_file do |dir|
      Misc.in_dir dir do
        CMD.cmd("wget http://www.thua45.cn/geredb/GREDB_dump.rar -O file.rar")
        CMD.cmd("unrar x file.rar")
        tsv = TSV.open 'Links.txt', :header_hash => '', :type => :double, :merge => true
        evidence = TSV.open "Evidence.txt", :header_hash => '', :type => :double, :key_field => 'lid', :merge => true
        evidence = evidence.to_list{|l| l* ";"}
        tsv.attach evidence, :merge => true, :fields => %w(pmid stn)
      end
    end
    tsv.to_s
  end

  GEREDB.claim GEREDB.tf_tg, :proc do 
    tsv = GEREDB[".source"].data.tsv :key_field => "reg", :fields => %w(tar eff pmid stn), :type => :double, :namespace => GEREDB.organism, :merge => true
    tsv.key_field = "Transcription Factor (Associated Gene Name)"
    tsv.fields = ["Target Gene (Associated Gene Name)", "Effect", "PMID", "Sentence Number"]

    dbtfs = TFClass.tfs.list

    tsv.select :key => dbtfs
  end


end

iif GEREDB.tf_tg.produce(true).find if __FILE__ == $0

