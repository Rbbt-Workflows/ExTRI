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

  TFClass.claim TFClass.hierarchy, :proc do
    url = "http://tfclass.bioinf.med.uni-goettingen.de/suppl/tfclass.ttl.gz"
		code_labels = {}
		Open.open(url) do |f|
			while line = f.gets
				line = Misc.fixutf8  line
				if m = line.match(%r(http://sybig.de/tfclass#((?:(?:\d+)\.)*(?:\d+))> rdf))
					codes = m.captures[0].split(".")
					l1 = f.gets
					l2 = f.gets
					l3 = f.gets
					l2 = Misc.fixutf8  l2
					label = l2.match(/"(.*)"/).captures[0]

          code_labels[codes*"."] = label
				end
			end
		end
    tsv = TSV.setup(code_labels, :key_field => "Term", :fields => ["Label"], :type => :single)
    tsv.to_s
  end

  def self.build_hier(hier, code)
    subset = "Class"
    info = {:name => hier[code][:name], :id => code, :subset => subset}
    if hier[code][:children]
      children = hier[code][:children].collect{|c|
        build_hier(hier, c)
      }
      info[:children] = children
    end
    if code.split(".").length == 5
      info[:size] = 1
    end
    info
  end

  TFClass.claim TFClass.hierarchy_json, :proc do 
    tsv = TFClass.hierarchy.tsv

    hier = {}
    first = []
    tsv.each do |term, label|
      codes = term.split(".")
      if codes.length == 1
        first << codes * "."
      end
      pre = codes[0..-2]
      hier[codes*"."] ||= {}
      hier[codes*"."][:name] = label
      hier[pre*"."] ||= {}
      hier[pre*"."][:children] ||= []
      hier[pre*"."][:children] << codes * "."
    end
    hier[''] = {:name => 'DbTFs', :children => first}
    
    build_hier(hier, '').to_json

  end

  TFClass.claim TFClass.tf_genus, :proc do 
    uni_equivalences = PRO.uniprot_equivalences.tsv :merge => true, :persist => true, :type => :flat
    uni2name = Organism.identifiers(ExTRI.organism).index :target => "Associated Gene Name", :persist => true
    gene2uniHsa = Organism.identifiers(Intact.organism("Hsa")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniMmu = Organism.identifiers(Intact.organism("Mmu")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true
    gene2uniRno = Organism.identifiers(Intact.organism("Rno")).index :target => "UniProt/SwissProt Accession", :order => true, :persist => true

    code_genus = {}
    Rbbt.share.databases.ExTRI.Nov2017_update["TFClass"].glob("*.tsv").each do |file|
      TSV.traverse file, :type => :array do |line|
        genus, code = line.split("\t")
        genus.sub!('"','')
        codes = [code, gene2uniHsa[code], gene2uniMmu[code], gene2uniRno[code]].compact.uniq
        codes += codes.collect{|code| uni_equivalences[code]}.compact.flatten.uniq
        names = uni2name.values_at(*codes.uniq).compact

        code_genus[names.first] = genus
        all_names << names.first unless names.empty?
      end
    end

    TSV.traverse Rbbt.share.databases.ExTRI.Nov2017_update["TFClass"]["tfclasscode2genesymbol_dic.txt"], :type => :array do |line|
        genus, code = line.split(":")
        genus.sub!('"','')
        name = code.gsub('"','').gsub(',','')
        next unless uni2name.include? name
        code_genus[name] = genus
    end

    tsv = TSV.setup(code_genus, :key_field => "TF", :fields => ["Genus"], :type => :single)
    tsv.to_s
  end
  
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
  iif TFClass.hierarchy_json.produce(true).find 
end

