
Workflow.require_workflow "Graph"
Workflow.require_workflow "Genomics"
add_workflow Graph, true

Genomics.knowledge_base.register :TF, Rbbt.data.TF, :source => "TF Associated Gene Name=~Associated Gene Name", :target => "TG Associated Gene Name=~Associated Gene Name", :fields => ["Sentence", "Sentence score", "Interaction score"]


require 'helpers/graph'
require 'rbbt/rest/web_tool'
include Sinatra::RbbtToolHelper
