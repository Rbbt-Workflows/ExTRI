require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/sources/organism'
require 'rbbt/sources/PRO'

module TFClass
  extend Resource
  self.subdir = 'share/databases/TFClass'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  #TFClass.claim TFClass.tf_tg, :proc do 
  #  url="http://tfclass.bioinf.med.uni-goettingen.de/suplementary/TFClass_ontologies.zip"
  #end
  
  TFClass.claim TFClass.tfs, :proc do 
    uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
    uni2name = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :persist => true
    gene2uniHsa = Organism.identifiers(Intact.organism("Hsa")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniMmu = Organism.identifiers(Intact.organism("Mmu")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniRno = Organism.identifiers(Intact.organism("Rno")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true

    all_names = Set.new
    Rbbt.share.databases.ExTRI.Nov2017_update["TFClass"].glob("*.tsv").each do |file|
      TSV.traverse file, :type => :array do |line|
        genus, code = line.split("\t")
        codes = [code, gene2uniHsa[code], gene2uniMmu[code], gene2uniRno[code]].compact.uniq
        codes += codes.collect{|code| uni_equivalences[code]}.compact.flatten.uniq
        names = uni2name.values_at(*codes.uniq).compact

        all_names << names.first unless names.empty?
      end
    end

    TSV.traverse Rbbt.share.databases.ExTRI.Nov2017_update["TFClass"]["tfclasscode2genesymbol_dic.txt"], :type => :array do |line|
        genus, code = line.split(":")
        name = code.gsub('"','').gsub(',','')
        next unless uni2name.include? name
        all_names << name
    end

    all_names.to_a.sort * "\n" + "\n"
  end
end
if __FILE__ == $0
  require 'rbbt/workflow'
  Workflow.require_workflow "ExTRI"
  iif TFClass.tfs.produce(true).find 
end

