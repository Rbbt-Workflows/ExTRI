require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/FNL'
require 'rbbt/sources/HTRI'
require 'rbbt/sources/TRRUST'
require 'rbbt/sources/tfacts'
require 'rbbt/sources/TFCheckpoint'

module FNL
  extend Workflow
end

require 'rbbt/tasks/FNL/clean'
require 'rbbt/tasks/FNL/validation'
require 'rbbt/tasks/FNL/database_coverage'
require 'rbbt/tasks/FNL/ner'
require 'rbbt/tasks/FNL/statistics'
require 'rbbt/tasks/FNL/year'
require 'rbbt/tasks/FNL/psicquic'
require 'rbbt/tasks/FNL/biogateway'

if __FILE__ == $0
  require 'rbbt/util/R'
  ENV["DISPLAY"] = "localhost:12.0"

  Log.severity = 0
  tsv = FNL.job(:sentence_coverage).run
  present_fields = tsv.fields.select{|f| f =~ /present/}.reject{|f| f =~ /Encode/ }
  raise "STOP"

  analysis = 6
  case analysis
  when 1
    scores = tsv.column("Interaction score").to_list
    scores.fields = ["score"]

    scores.add_field "present" do |k,v|
      tsv[k].values_at(*present_fields).reject{|v| v.empty?}.any? ? "present" : "not present"
    end

    scores.R_interactive <<-EOF

    data = rbbt.tsv(data_file)
    data$present = as.factor(data$present)

    EOF
  when 1.1
    tsv = FNL.job(:sentence_coverage).run
    scores = tsv.column("Sentence score").to_list
    scores.fields = ["score"]

    scores.add_field "present" do |k,v|
      tsv[k].values_at(*present_fields).reject{|v| v.empty?}.any? ? "present" : "not present"
    end

    scores.R_interactive <<-EOF

    data = rbbt.tsv(data_file)
    data$present = as.factor(data$present)

    EOF
  when 1.2
    tsv = FNL.job(:sentence_coverage).run
    scores = tsv.slice(["Sentence score", "Interaction score"])
    Log.tsv scores
    scores.fields = ["score_sent", "score_int"]

    pmid_counts = TSV.setup({}, :key_field => 'Pair', :fields => ["Counts"], :type => :single)
    pmid_present = TSV.setup({}, :key_field => 'Pair', :fields => ["Present"], :type => :single)

    tsv.through do |k, values|
      tf, tg = values.values_at("TF Associated Gene Name", "TG Associated Gene Name")
      pair = [tf, tg] * ":"
      present = values.values_at(*present_fields).reject{|v| v.empty?}.any?
      pmid_counts[pair] ||= 0
      pmid_counts[pair] += 1
      pmid_present[pair] ||= true if present
    end

    scores.add_field "pmid_counts" do |k,v|
      pair = k.split(":").values_at(2,3) * ":"
      pmid_counts[pair] || 0
    end

    scores.add_field "pmid_present" do |k,v|
      pair = k.split(":").values_at(2,3) * ":"
      pmid_present[pair] || false
    end

    scores.R_interactive <<-EOF

    data = rbbt.tsv(data_file)
    model = lm(score_sent ~ score_int, data = data)
    r = summary(model)$r.squared

    plot(score_sent ~ score_int, data = data, main=paste("Sentence vs. Interaction scores. R-squared: ", r))
    lines(abline(model))

    EOF
  when 2
    pmid_counts = TSV.setup({}, :key_field => 'Pair', :fields => ["Counts"], :type => :single)
    pmid_present = TSV.setup({}, :key_field => 'Pair', :fields => ["Present"], :type => :single)

    tsv.through do |k, values|
      tf, tg = values
      pair = [tf, tg] * ":"
      present = values.values_at(*present_fields).reject{|v| v.empty?}.any?
      pmid_counts[pair] ||= 0
      pmid_counts[pair] += 1
      pmid_present[pair] ||= true if present
    end

    pmid_counts = pmid_counts.to_list
    pmid_counts.attach pmid_present.to_list

    pmid_counts.R_interactive <<-EOF

    data = rbbt.tsv(data_file)
    data$present = as.factor(data$present)

    EOF
  when 3
    triplet_counts = TSV.setup({}, :key_field => 'Triplet', :fields => ["Counts"], :type => :single)
    triplet_present = TSV.setup({}, :key_field => 'Triplet', :fields => ["Present"], :type => :single)

    tsv.through do |k, values|
      tf, tg = values
      pmid = k.split(":").first
      triplet = [pmid,tf, tg] * ":"
      present = values.values_at(*present_fields).reject{|v| v.empty?}.any?
      triplet_counts[triplet] ||= 0
      triplet_counts[triplet] += 1
      triplet_present[triplet] ||= true if present
    end

    triplet_counts = triplet_counts.to_list
    triplet_counts.attach triplet_present.to_list

    triplet_counts.R_interactive <<-EOF

    data = rbbt.tsv(data_file)
    data$present = as.factor(data$present)

    EOF
  when 4

    tsv = FNL.job(:all).run
    good_present_fields = tsv.fields.select{|f| f =~ /present/}.reject{|f| f =~ /Encode|GOA|TFChe/ }
    pair_present = tsv.slice(good_present_fields)


    pair_present.R_interactive <<-EOF

    data = rbbt.tsv(data_file)
    data = apply(data,2,function(x){ x != ""})
    fields = #{R.ruby2R good_present_fields}
    categories = #{R.ruby2R good_present_fields.collect{|f| num = pair_present.select(f){|l| ! (l.nil? || l.empty?)}.length; [f, num] * ": "}}
    rbbt.plot.venn(data, fields, category=categories, margin=0.2)

    EOF
  when 5

    tsv = FNL.job(:all).path.tsv :key_field => "TF", :merge => true, :zipped => true
    good_present_fields = tsv.fields.select{|f| f =~ /present/}.reject{|f| f =~ /Encode|GOA|TFChe/ }
    tf_present = tsv.slice(good_present_fields).to_list{|v| v.compact.uniq.first}

    tf_present.R_interactive <<-EOF

    data = rbbt.tsv(data_file)
    data = apply(data,2,function(x){ x != ""})
    fields = #{R.ruby2R good_present_fields}
    categories = #{R.ruby2R good_present_fields.collect{|f| num = tf_present.select(f){|l| ! (l.nil? || l.empty?)}.length; [f, num] * ": "}}
    rbbt.plot.venn(data, fields, category=categories, margin=0.2)

    EOF

  when 6

    validation = Rbbt.data["Astrid Validation"]["validation_samples_260317.tsv"].tsv :type => :list
    validation = validation.attach tsv, :fields => ["Interaction score"]

    all = FNL.job(:all).run
    good_present_fields = all.fields.select{|f| f =~ /present/}.reject{|f| f =~ /Encode|GOA|TFChe/ }
    pair_present = all.slice(good_present_fields)
    pair_present.fields = pair_present.fields.collect{|f| name = f.match(/\[(.*)\]/)[1]; "pair #{name}"}

    validation.add_field pair_present.key_field do |k,v|
      [v[0], v[1]] * ":"
    end
    validation.attach pair_present

    all_tf = FNL.job(:all).path.tsv :key_field => "TF", :merge => true, :zipped => true
    good_present_fields = all_tf.fields.select{|f| f =~ /present/}.reject{|f| f =~ /Encode|GOA|TFChe/ }
    tf_present = all_tf.slice(good_present_fields).to_list{|v| v.compact.uniq.first}
    tf_present.fields = tf_present.fields.collect{|f| name = f.match(/\[(.*)\]/)[1]; "TF #{name}"}
    validation.add_field tf_present.key_field do |k,v|
      v[0]
    end
    validation.attach tf_present

    validation.add_field "Triplet" do |k,v|
      k.split(":").values_at(0,2,3) * ":"
    end

    validation.add_field "Valid" do |k,v|
      v["1"] and not v["1"] == ""
    end

    validation.add_field "Valid ext" do |k,v|
      (v["1"] and not v["1"] == "") or (v["0"] == "n")
    end

    pmid_counts = TSV.setup({}, :key_field => 'TF:TG', :fields => ["PMID Counts"], :type => :single)
    pmid_present = TSV.setup({}, :key_field => 'TF:TG', :fields => ["Present"], :type => :single)

    sentence_counts = TSV.setup({}, :key_field => 'Triplet', :fields => ["Sentence Counts"], :type => :single)
    tsv.through do |k, values|
      tf, tg = values
      pair = [tf, tg] * ":"
      triplet = k.split(":").values_at(0,2,3) * ":"
      present = values.values_at(*present_fields).reject{|v| v.empty?}.any?
      pmid_counts[pair] ||= 0
      pmid_counts[pair] += 1
      pmid_present[pair] ||= true if present
      pmid_present[pair] ||= true if present
      sentence_counts[triplet] ||= 0
      sentence_counts[triplet] += 1
    end

    validation.attach pmid_counts.to_list
    validation.attach pmid_present.to_list
    validation.attach sentence_counts.to_list

    sent = FNL.job(:sentence_coverage).run
    validation.attach sent, :fields => "Sentence score"

    Open.write('/tmp/sentence_validation.tsv', validation.to_s)
    raise "STOP"

    validation.R_interactive <<-EOF
data = rbbt.tsv(data_file)
names(data) <- make.names(names(data))
    EOF

  end


end

#require 'MODULE/tasks/basic.rb'

#require 'rbbt/knowledge_base/FNL'
#require 'rbbt/entity/MODULE'

