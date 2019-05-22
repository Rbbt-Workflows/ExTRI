require 'rbbt/tsv/excel'
module ExTRI

  dep :ExTRI_confidence
  extension :tsv
  task :info_for_training => :text do
    conf = step(:ExTRI_confidence).load.keys 
    tsv = TSV.excel(Rbbt.data.training_subset["pairs_261118_for_training.xlsx"])
    
    dumper = TSV::Dumper.new :key_field => "PMID", :fields => ["Sentence ID", "Pair", "Code"] + tsv.fields, :type => :list
    dumper.init

    TSV.traverse tsv, :bar => true, :into => dumper do |pair, values|
      matches = conf.select{|k| k.include? ":" + pair}
      res = matches.collect do |m|
        pmid,num,tf,tg = m.split(":")
        [pmid, [num, pair, m] + values]
      end

      res.extend MultipleResult
      res
    end
  end

  dep :ExTRI_confidence
  task :nuclear_receptor_count => :tsv do
    receptors = TFClass.hierarchy.produce.read.split("\n").select{|l| l =~ /^2\.1\.\d\.\d\.\d/}.collect{|l| l.split("\t").last.split(/[\s(), ]+/)}.flatten
    counts = TSV.setup({},"TF~Count#:type=:single")
    tsv = step(:ExTRI_confidence).load
    TSV.traverse tsv.select("Transcription Factor (Associated Gene Name)" => receptors) do |key,values|
      tf,tg, *rest = values
      counts[tf] ||= 0
      counts[tf] += 1
    end
    counts
  end

  dep :sentence_coverage
  task :DB_signed_overlap => :tsv do
    tsv = step(:sentence_coverage).load
    overlap = tsv.select("[TRRUST] Regulation"){|v| (v.split(";") & %w(Activation Repression)).any?}
    overlap.merge!(tsv.select("[TFacts] Sign"){|v| (v.split(";") & %w(UP DOWN)).any?})
    overlap
  end

  dep :DB_signed_overlap
  task :nuclear_factor_signed_overlap => :tsv do
    receptors = TFClass.hierarchy.produce.read.split("\n").select{|l| l =~ /^2\.1\.\d\.\d\.\d/}.collect{|l| l.split("\t").last.split(/[\s(), ]+/)}.flatten
    tsv = step(:DB_signed_overlap).load
    tsv.select("Transcription Factor (Associated Gene Name)" => receptors)
  end

  dep :DB_signed_overlap
  extension :tsv
  task :not_nuclear_factor_signed_overlap => :tsv do
    receptors = TFClass.hierarchy.produce.read.split("\n").select{|l| l =~ /^2\.1\.\d\.\d\.\d/}.collect{|l| l.split("\t").last.split(/[\s(), ]+/)}.flatten
    tsv = step(:DB_signed_overlap).load
    tsv = tsv.select({"Transcription Factor (Associated Gene Name)" => receptors}, true)
    fabios = Rbbt.sandbox.sentence_id.tsv.values.flatten.compact
    iif fabios
    tsv.select(:key){|k| ! fabios.include? k.split(":").values_at(0,1) * ":"}
  end
end
