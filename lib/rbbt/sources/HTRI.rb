require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'

module HTRI
  extend Resource
  self.subdir = 'share/databases/HTRI'

  def self.organism(org="Hsa")
    Organism.default_code(org)
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  HTRI.claim HTRI[".source"].data, :proc do |filename|
    raise "Please download http://www.lbbc.ibb.unesp.br/htri/datasHTRIdb?down=1 from http://www.lbbc.ibb.unesp.br/htri/index.jsp and place it in #{filename}"
  end

  HTRI.claim HTRI.tf_tg, :proc do 
    file = HTRI['.source'].data
    tsv = TSV.open(file, :fix => Proc.new{|l| l.gsub(/\t+/,"\t")}, :header_hash => "", :key_field => "SYMBOL_TF", :fields => ["SYMBOL_TG", "TECHNIQUE", "PUBMED_ID"], :merge => true)
    tsv.key_field = "Transcription Factor (Associated Gene Name)"
    tsv.fields = ["Target Gene (Associated Gene Name)", "Technique", "PMID"]
    tsv.namespace = HTRI.organism

    tsv.add_field "Confidence" do |k,v|
      v["TECHNIQUE"].collect{|technique|
        case technique
        when "Electrophoretic Mobility Shift Assay", "Chromatin Immunoprecipitation", "Dnase I Footprinting", 
          "Avidin/Biotin-conjugated DNA Binding Assay", "CpG Chromatin Immunoprecipitation", "Surface Plasmon Resonance",
          "Concatenate Chromatin Immunoprecipitation", "DNA Affinity Precipitation Assay", 
          "Yeast One-Hybrid Assay", "Southwestern Blotting", "DNA Affinity Chromatography", 
          "Streptavidin Chromatin Immunoprecipitation" 

          "High"

        when "Chromatin Immunoprecipitation coupled with microarray", "Chromatin Immunoprecipitation coupled with deep sequencing"
          "Low"
        else
          Log.warn "Unknown technique in HTRI: " << technique
          "Unknown"
        end
      }
    end

    tsv.to_s
  end
end

Log.tsv HTRI.tf_tg.produce(true).tsv if __FILE__ == $0


