#!/usr/bin/env ruby

require 'rbbt-util'
require 'rbbt/util/simpleopt'

$0 = "rbbt #{$previous_commands*""} #{ File.basename(__FILE__) }" if $previous_commands

options = SOPT.setup <<EOF

Matche HTRIDB to TRI tagged pairs

Add column matching small_scale_experiment_only.csv from HTRIDB to the
sentences.TRI_tagged.pairs_cancer_low_genes_disease.tsv file

$ #{$0} <sentences.TRI_tagged.pairs_cancer_low_genes_disease.tsv> <small_scale_experiment_only.csv>

-h--help Print this help

EOF

error = true unless ARGV.length == 2

if options[:help] or error
  if defined? rbbt_usage
    rbbt_usage 
  else
    puts SOPT.doc
  end
  exit 0
end

tri, sse = ARGV

sse = TSV.open sse, :sep => ';', :header_hash => '', :type => :list, :fields => ["SYMBOL_TF","SYMBOL_TG"]

pairs = Set.new
TSV.traverse sse, :into => pairs do |code,values|
  tf,tg = values
  [tf,tg] * "~"
end

key_field = "Sentence ID"
fields = ["TF", "TG", "Score", "Sentence", "Type", "Match key", "Interaction key", "Other", "HTRIDB"]


dumper = TSV::Dumper.new :key_field => key_field, :fields => fields, :type => :list, :merge => true
dumper.init
TSV.traverse tri, :type => :list, :into => dumper do |id,values|
  id = id.first if Array === id
  tf,tg,score,sentence,type,match,key,other = values

  pair = [tf,tg] * "~"
  c = pairs.include?(pair) ? "YES" : "NO"
  new_values = [tf,tg,score,sentence,type,match,key,other, c]
  [id,new_values]
end

Misc.consume_stream dumper.stream, false, STDOUT
