require 'rbbt/knowledge_base'
require 'rbbt/sources/organism'

module FNL

  class << self 
    attr_accessor :knowledge_base_dir
  end
  self.knowledge_base_dir = Rbbt.var.knowledge_base.FNL

  def self.organism
    Organism.default_code("Hsa")
  end

  def self.knowledge_base
    @knowledge_base ||= begin
                          kb = KnowledgeBase.new self.knowledge_base_dir, self.organism

                          kb
                        end
  end
end

FNL.knowledge_base.register :TF, Rbbt.data.TF, :source => "TF Associated Gene Name=~Associated Gene Name", :target => "TG Associated Gene Name=~Associated Gene Name", :fields => ["Sentence", "Sentence score", "Interaction score"]
