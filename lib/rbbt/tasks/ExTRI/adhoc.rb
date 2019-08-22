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
  
  dep :pairs
  task :cyt_reg_analysis => :tsv do
    pairs = step(:pairs).load
    cyt_reg = pairs.select("[CytReg] present" => "CytReg")
    trrust = pairs.select("[TRRUST] present" => "TRRUST")

    cyt_reg_pmids = cyt_reg.column("[CytReg] PMIDs").values.flatten.collect{|v| v.split(";")}.flatten.uniq
    trrust_pmids = trrust.column("[TRRUST] PMID").values.flatten.collect{|v| v.split(";")}.flatten.uniq

    pmids = cyt_reg_pmids & trrust_pmids

    pmids_pairs = TSV.setup({}, "PMID~CytReg,TRRUST,Common,CytReg Signed,TRRUST Signed,Common Signed#:type=:double")
    pmids.each do |pmid|
      pairs_cr = cyt_reg.select("[CytReg] PMIDs" => /\b#{pmid}\b/).keys
      pairs_t = trrust.select("[TRRUST] PMID" => /\b#{pmid}\b/).keys

      pairs_cr_s = pairs_cr.collect{|p| cyt_reg[p]["[CytReg] Activation/Repression"].collect{|a| [p,a]*"~"}}.flatten.uniq
      pairs_t_s = pairs_t.collect{|p| trrust[p]["[TRRUST] Regulation"].collect{|a| [p,a]*"~"}}.flatten.uniq
      pairs_t_s.reject!{|p| p.include? "Unknown"}

      common = pairs_cr & pairs_t
      common_s = pairs_cr_s & pairs_t_s
      pmids_pairs[pmid] = [pairs_cr, pairs_t, common, pairs_cr_s, pairs_t_s, common_s]
    end

    pmids_pairs
  end

  dep :cyt_reg_analysis
  task :cyt_reg_statistics => :tsv do

    tsv = step(:cyt_reg_analysis).load

    res = TSV.setup({}, "PMID~All,CytReg,TRRUST,Common,Overlap %,Found by CytReg %,Found by TRRUST %,All Signed,CytReg Signed,TRRUST Signed,Common Signed,Overlap % Signed,Found by CytReg % Signed,Found by TRRUST Signed %#:type=:list#:cast=:to_f")

    tsv.each do |pmid,values|
      cr, t, c, cr_s, t_s, c_s = values.collect{|v| v.length}

      all = cr + t - c
      overlap = c.to_f / all
      found_by_cr = cr.to_f / all
      found_by_t = t.to_f / all

      all_s = cr_s + t_s - c_s
      overlap_s = c_s.to_f / all_s
      found_by_cr_s = cr_s.to_f / all_s
      found_by_t_s = t_s.to_f / all_s


      res[pmid] = [all, cr, t, c, (overlap * 100).to_i, (found_by_cr * 100).to_i, (found_by_t * 100).to_i, all_s, cr_s, t_s, c_s, (overlap_s * 100).to_i, (found_by_cr_s * 100).to_i, (found_by_t_s * 100).to_i]
    end
    
    res
  end
end
