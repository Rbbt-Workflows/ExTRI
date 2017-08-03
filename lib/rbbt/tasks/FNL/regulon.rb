module FNL

  dep :pairs
  input :tfs, :array, "TFs to consider"
  input :consensus, :integer, "Minumum number of DBs containing pair", 2
  input :high, :boolean, "Consider only high-confidence", true
  input :sample, :float, "How many targets to select randomly", 1
  task :regulon => :tsv do |tfs,consensus,high,sample|

    dumper = TSV::Dumper.new(:key_field => "TF (Associated Gene Name)", :fields => ["TG (Associated Gene Name)", "Sign"], :type => :double)
    dumper.init

    parser = TSV::Parser.new step(:pairs)
    fields = parser.fields
    TSV.traverse parser, :into => dumper do |pair, values|
      tf, tg = pair.split(":")
      next unless tfs.nil? or tfs.include? tf
      skip = []
      present = []
      not_tf_class = []
      signs = []
      fields.zip(values).each do |field,value|
        skip << field if high and field.include?("confidence") and value == "Low"
        present << field if field.include?("present") and value != ""
        not_tf_class << field if field.include?("TFClass_human") and value != "TFclass"
        signs += value.split(";") if field.include?("Sign")
        signs += value.split(";") if field.include?("Regulation")
      end
      present = present.collect{|f| f.match(/\[(.*)\]/)[1]}
      skip = skip.collect{|f| f.match(/\[(.*)\]/)[1]}

      if signs.include?("Repression") or signs.include?("DOWN")
        sign = -1
      elsif signs.include?("Activation") or signs.include?("UP")
        sign = 1
      else
        sign = nil
      end

      next if not_tf_class.any? and not present.include? "TFacts"
      next unless (present - skip).length >= consensus
      [tf, [tg, sign]]
    end
  end
end
