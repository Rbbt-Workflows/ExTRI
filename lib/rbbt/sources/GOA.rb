require 'rbbt/sources/go'
require 'rbbt/sources/PRO'

module GO
  #GO.claim GO.tf_tg, :proc do 

  #  tsv = TSV.setup({}, :key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "Sign", "GO Term"], :type => :double, :namespace => FNL.organism)
  #  
  #  mgi2uni = Organism.identifiers("Mmu/feb2014").index :target => "UniProt/SwissProt Accession", :persist => true

  #  uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
  #  uni2name = Organism.identifiers(FNL.organism).index :target => "Associated Gene Name", :persist => true

  #  parser = TSV::Parser.new Rbbt.share.databases.FNL.Nov2017_update['GOA.Mmu-Hsa.tsv'], :type => :list, :header_hash => ''
  #  fields = parser.fields
  #  TSV.traverse parser, :into => tsv, :bar => true do |key, parts|
  #    parts = NamedArray.setup(parts, fields)

  #    res = [].extend MultipleResult

  #    tf_uni = parts["GENE PRODUCT ID"]
  #    with = parts["WITH/FROM"]
  #    term = parts["GO TERM"]
  #    name = parts["GO NAME"]
  #    ext = parts["ANNOTATION EXTENSION"]

  #    tgs_uni = ext.split("|").select{|f| f.include? "UniProtKB"}.collect{|f| f.split(":").last[0..-2]}
  #    tgs_mgi = ext.split("|").select{|f| f.include? "MGI"}.collect{|f| "MGI:" + f.split(":").last[0..-2]}

  #    tgs_uni += mgi2uni.values_at(*tgs_mgi).compact

  #    next if tgs_uni.empty?

  #    tf_ext_uni = uni_equivalences[tf_uni] 
  #    tgs_ext_uni = uni_equivalences.values_at(*tgs_uni)

  #    tf = uni2name.values_at(*tf_ext_uni).compact.first
  #    tgs_ext_uni.each do |tg_ext_uni|
  #      tg = uni2name.values_at(*tg_ext_uni).compact.first

  #      next if tg.nil? or tf.nil?

  #      sign = "DOWN" if name.include?("negative regulation") or name.include?("transcriptional repressor activity")
  #      sign = "UP" if name.include?("positive regulation") or name.include?("transcriptional activator activity")
  #      sign = "Unknown" if sign.nil?

  #      res << [tf, [tg, sign, term]]
  #    end

  #    res.uniq!

  #    res
  #  end

  #  tsv
  #end
  GO.claim GO.tf_tg, :proc do 
    tsv = TSV.setup({}, :key_field => "Transcription Factor (Associated Gene Name)", :fields => ["Target Gene (Associated Gene Name)", "Sign"], :type => :double, :namespace => FNL.organism)

    uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
    uni2name = Organism.identifiers(FNL.organism).index :target => "Associated Gene Name", :persist => true
    gene2uniHsa = Organism.identifiers(Intact.organism("Hsa")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniMmu = Organism.identifiers(Intact.organism("Mmu")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniRno = Organism.identifiers(Intact.organism("Rno")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true

    Rbbt.share.databases.FNL.Nov2017_update["GOA.source"].glob("*.tsv").each do |file|
      TSV.traverse file, :type => :array do |line|
        tf_uni, tf_name, sign, relations_str = line.split("\t")

        tf_unis = [tf_uni] + (uni_equivalences[tf_uni] || [])
        tf_name = uni2name.values_at(*tf_unis).compact.first

        tg_names = []
        relations_str.split(/[,|]/).select do |rel|
          type, gene_str = rel.split("(")
          next unless %w(has_direct_input has_input has_regulation_target regulates_transcription_of).include? type
          code = gene_str.partition(":").last[0..-2]
          codes = [code, gene2uniHsa[code], gene2uniMmu[code], gene2uniRno[code]].compact.uniq
          codes += codes.collect{|code| uni_equivalences[code]}.compact.flatten.uniq
          names = uni2name.values_at(*codes.uniq).compact
          tg_names << names.first
        end


        sign = nil if sign == "NA"
        next if tf_name.nil?
        tg_names.compact.each do |tg_name|
          tsv.zip_new tf_name, [tg_name, sign]
        end
      end
    end

    tsv.remove_duplicates.to_s
  end
end

if __FILE__ == $0
  require 'rbbt/workflow'
  Workflow.require_workflow "FNL"
  iii GO.tf_tg.produce(true).find
end
