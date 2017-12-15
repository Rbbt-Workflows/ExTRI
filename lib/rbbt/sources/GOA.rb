require 'rbbt/sources/go'
require 'rbbt/sources/PRO'

module GO
  GO.claim GO.tf_tg, :proc do 

    tsv = TSV.setup({}, :key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "Sign", "GO Term"], :type => :double, :namespace => FNL.organism)
    
    mgi2uni = Organism.identifiers("Mmu/feb2014").index :target => "UniProt/SwissProt Accession", :persist => true

    uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
    uni2name = Organism.identifiers(FNL.organism).index :target => "Associated Gene Name", :persist => true

    parser = TSV::Parser.new Rbbt.share.databases.FNL.Nov2017_update['GOA.Mmu-Hsa.tsv'], :type => :list, :header_hash => ''
    fields = parser.fields
    TSV.traverse parser, :into => tsv, :bar => true do |key, parts|
      parts = NamedArray.setup(parts, fields)

      res = [].extend MultipleResult

      tf_uni = parts["GENE PRODUCT ID"]
      with = parts["WITH/FROM"]
      term = parts["GO TERM"]
      name = parts["GO NAME"]
      ext = parts["ANNOTATION EXTENSION"]

      tgs_uni = ext.split("|").select{|f| f.include? "UniProtKB"}.collect{|f| f.split(":").last[0..-2]}
      tgs_mgi = ext.split("|").select{|f| f.include? "MGI"}.collect{|f| "MGI:" + f.split(":").last[0..-2]}

      tgs_uni += mgi2uni.values_at(*tgs_mgi).compact

      next if tgs_uni.empty?

      tf_ext_uni = uni_equivalences[tf_uni] 
      tgs_ext_uni = uni_equivalences.values_at(*tgs_uni)

      tf = uni2name.values_at(*tf_ext_uni).compact.first
      tgs_ext_uni.each do |tg_ext_uni|
        tg = uni2name.values_at(*tg_ext_uni).compact.first

        next if tg.nil? or tf.nil?

        sign = "DOWN" if name.include?("negative regulation") or name.include?("transcriptional repressor activity")
        sign = "UP" if name.include?("positive regulation") or name.include?("transcriptional activator activity")
        sign = "Unknown" if sign.nil?

        res << [tf, [tg, sign, term]]
      end

      res.uniq!

      res
    end

    tsv
  end
end

if __FILE__ == $0
  require 'rbbt/workflow'
  Workflow.require_workflow "FNL"
  iii GO.tf_tg.produce(true).find
end
