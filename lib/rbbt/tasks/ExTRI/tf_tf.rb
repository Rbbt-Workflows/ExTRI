module ExTRI

  dep :ExTRI_confidence, :pmids => 2, :sentences => 2, :score => 1.6, :test_set => []
  task :tf_tf => :tsv do
    tsv = step(:ExTRI_confidence).load
    tfs = tsv.column("Transcription Factor (Associated Gene Name)").values
    tfs -= %w(AP1 NFKB)
    tsv = tsv.select("Target Gene (Associated Gene Name)" => tfs)
    tsv = tsv.select("Transcription Factor (Associated Gene Name)" => tfs)

    tsv.add_field "Auto-regulation" do |p,values|
      values[1] == values[0] ? "Auto" : ""
    end
  end
 
  dep :tf_tf
  task :tf_tf_counts => :tsv do
    tsv = TSV.setup({}, "TF:TG~Triplets#:type=:flat")
    TSV.traverse step(:tf_tf) do |triplet, values|
      tf, tg, *rest = values
      pair = [tf,tg].sort * ":"
      if tsv[pair].nil?
        tsv[pair] = [triplet]
      else
        tsv[pair] << triplet
      end
    end

    tsv = tsv.to_double

    tsv.add_field "TF" do |pair, values|
      pair.split(":").first
    end

    tsv.add_field "TG" do |pair, values|
      pair.split(":").last
    end

    tsv.add_field "Auto-regulation" do |pair, values|
      pair.split(":").uniq.length == 1 ? "Auto" : ""
    end

    tsv.add_field "PMID (count)" do |pair, values|
      pmids = values.first.collect{|trip| trip.split(":").first}
      pmids.uniq.length
    end

    tsv.add_field "Sentences (count)" do |pair, values|
      sentences = values.first.collect{|trip| trip.split(":").values_at(0,1) * ":"}
      sentences.uniq.length
    end

   tsv.reorder :key, ["TF", "TG", "Auto-regulation", "PMID (count)", "Sentences (count)", "Triplets"]
  end

  dep :ExTRI_confidence, :pmids => 2, :sentences => 2, :score => 1.6, :test_set => []
  input :same_tg, :boolean, "Consider only co-tfs for the same tf", false
  task :co_tf => :tsv do |same_tg|
    fnl = step(:ExTRI_confidence).load
    article_codes = TSV.setup({}, :key_field => "PMID", :fields => ["Codes"], :type => :flat)
    sentence_codes = TSV.setup({}, :key_field => "Sentence", :fields => ["Codes"], :type => :flat)
    sentence_text = {}
    TSV.traverse fnl, :bar => true do |code, values|
      pmid, sentence_no, tf, tg = code.split(":")

      sentence = [pmid, sentence_no] * ":"
      article_codes[pmid] ||= []
      article_codes[pmid] << code
      sentence_codes[sentence] ||= []
      sentence_codes[sentence] << code
      sentence_text[sentence] ||= values[3]
    end

    article_co_tf = TSV.setup({}, :key_field => "PMID", :fields => ["Transcription Factor (Associated Gene Name)"], :type => :flat)
    TSV.traverse article_codes do |pmid, codes|
      tfs = codes.collect{|code| code.split(":")[2]}
      tfs -= %w(AP1 NFKB)
      next unless tfs.length > 1
      article_co_tf[pmid] = tfs.uniq
    end

    sentence_co_tf = TSV.setup({}, :key_field => "Sentence", :fields => ["Transcription Factor (Associated Gene Name)"], :type => :flat)
    TSV.traverse sentence_codes do |sentence, codes|
      tfs = codes.collect{|code| code.split(":")[2]}
      tfs -= %w(AP1 NFKB)
      next unless tfs.length > 1
      sentence_co_tf[sentence] = tfs.uniq
    end

    co_tf = TSV::Dumper.new(:key_field => "CO TF candidates", :fields => ["Type", "Code", "Text"], :type => :doble)
    co_tf.init

    Thread.new do 
      Misc.consume_stream co_tf.stream, false, tmp_path
    end

    article_co_tf.each do |pmid, codes|
      pairs = Set.new 
      codes.each do |c1|
        codes.each do |c2|
          next if c1 == c2
          pairs << ([c1,c2].sort * ":")
        end
      end

      pairs.to_a.each do |pair|
        co_tf.add pair, ['PMID', pmid, '']
      end
    end

    sentence_co_tf.each do |sentence, codes|
      pairs = Set.new 
      codes.each do |c1|
        codes.each do |c2|
          pairs << [c1,c2].sort * ":"
        end
      end

      pairs.to_a.each do |pair|
        co_tf.add pair, ['Sentence', sentence, sentence_text[sentence]]
      end
    end

    co_tf.close

    nil
  end

  dep :co_tf
  task :co_tf_counts => :tsv do
    tsv = TSV.setup({}, "Co-TF candidate", :fields => %w(PMID Sentence), :type => :list, :cast => :to_i)
    TSV.traverse step(:co_tf) do |pair, values|
      pair = pair.first
      type, code, text = values

      tsv[pair] ||= [0, 0]
      if type.first == "PMID"
        tsv[pair][0]+=1
      else
        tsv[pair][1]+=1
     end

    end
    tsv
  end
end
