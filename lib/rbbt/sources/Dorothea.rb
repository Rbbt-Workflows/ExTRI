require 'rbbt-util'
require 'rbbt/resource'

module Dorothea
  extend Resource
  self.subdir = 'share/databases/DoRothEA'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Dorothea.claim Dorothea.dorothea_a, :proc do 
    Rbbt.share.databases.ExTRI["DorotheaA.csv"].tsv :header_hash => "", :sep => ',', :merge => true

  end

  Dorothea.claim Dorothea.tf_tg, :proc do 
    dorotheaA = Dorothea.dorothea_a.tsv(:merge => true, :key_field => "source_genesymbol", :fields => %w(target_genesymbol references consensus_direction consensus_stimulation consensus_inhibition))
    dorotheaA.key_field = "Transcription Factor (Associated Gene Name)"
    dorotheaA.fields = ["Target Gene (Associated Gene Name)", "PMID", "Directed", "Stimulate", "Inhibit"]

    dorotheaA.process "Directed" do |v|
      v.collect{|e| e == "1" ? true : false}
    end

    dorotheaA.process "PMID" do |v|
      v.collect{|e| e == "NA" ? [] : e.split(";").collect{|v| v.split(":").last}.uniq * ";" }
    end

    dorotheaA.add_field "Effect" do |k,v|
      target, pmid, directed, stimulate, inhibit = v
      stimulate.zip(inhibit).collect{|p| 
        case p
        when %w(1 1)
          "Both"
        when %w(1 0)
          "Stimulate"
        when %w(0 1)
          "Inhibit"
        else
          "Unknown"
        end
      }
    end


    dorotheaA = dorotheaA.slice dorotheaA.fields - %w(Stimulate Inhibit)
  end
end

iif Dorothea.tf_tg.produce(true).find if __FILE__ == $0

