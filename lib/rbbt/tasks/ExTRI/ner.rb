module ExTRI

  dep :sentence_coverage
  task :sentence_coverage_NER => :tsv do

    ner = Rbbt.data.Martin_NER["ExTRI.sentence_coverage_full_genes.tsv.gz"].tsv :type => :double, :merge => true
    ner.key_field = "Sentence ID"
    ner.fields = ["?", "Text", "??", "Entrez Gene ID", "Tax ID", "Confidence"]
    ner = ner.slice(["Text", "Entrez Gene ID", "Tax ID", "Confidence"])

    name2entrez = Organism.identifiers(ExTRI.organism).index :target => "Entrez Gene ID", :order => true, :persist => true
    tsv = step(:sentence_coverage).load


    tsv.add_field "TF NER INFO" do |k,v|
      sentence_id = k.split(":").values_at(0,1) * ":"
      ner_values = ner[sentence_id]
      sentence = v["Sentence"]

      gene = v["Transcription Factor (Associated Gene Name)"]
      entrez = name2entrez[gene]
      if ner_values.nil? or entrez.nil?
        nil
      else
        matches = {}
        Misc.zip_fields(ner_values).each do |text, n_entrez, tax, conf|
          next if entrez != n_entrez
          start = sentence.index(text)
          eend = start + text.length 
          matches[tax] = [text,tax,start,eend] * "::"
        end
        matches["9606"] || matches.values.first
      end
    end

    tsv.add_field "TF Text" do |k,values|
      next if values["TF NER INFO"].nil?
      text, tax, start, eend = values["TF NER INFO"].split("::")
      text
    end

    tsv.add_field "TF Tax" do |k,values|
      next if values["TF NER INFO"].nil?
      text, tax, start, eend = values["TF NER INFO"].split("::")
      tax
    end

    tsv.add_field "TF range" do |k,values|
      next if values["TF NER INFO"].nil?
      text, tax, start, eend = values["TF NER INFO"].split("::")
      [start, eend] * ":"
    end

    tsv.add_field "TG NER INFO" do |k,v|
      sentence_id = k.split(":").values_at(0,1) * ":"
      ner_values = ner[sentence_id]
      sentence = v["Sentence"]

      gene = v["Target Gene (Associated Gene Name)"]
      entrez = name2entrez[gene]
      if ner_values.nil? or entrez.nil?
        nil
      else
        matches = {}
        Misc.zip_fields(ner_values).each do |text, n_entrez, tax, conf|
          next if entrez != n_entrez
          start = sentence.index(text)
          eend = start + text.length 
          matches[tax] = [text,tax,start.to_s,eend.to_s] * "::"
        end
        matches["9606"] || matches.values.first
      end
    end

    tsv.add_field "TG Text" do |k,values|
      next if values["TG NER INFO"].nil?
      text, tax, start, eend = values["TG NER INFO"].split("::")
      text
    end

    tsv.add_field "TG Tax" do |k,values|
      next if values["TG NER INFO"].nil?
      text, tax, start, eend = values["TG NER INFO"].split("::")
      tax
    end

    tsv.add_field "TG range" do |k,values|
      next if values["TG NER INFO"].nil?
      text, tax, start, eend = values["TG NER INFO"].split("::")
      [start, eend] * ":"
    end

    tsv.slice(tsv.fields - ["TF NER INFO", "TG NER INFO"])

  end

  #dep :sentence_coverage_full_NER
  #task :sentence_coverage_full_NER_counts => :tsv do
  #  tsv = step(:sentence_coverage_full_NER).load

  #  pmid_counts = TSV.setup({}, :key_field => 'Pair', :fields => ["Counts"], :type => :single)
  #  sentence_counts = TSV.setup({}, :key_field => 'Triplet', :fields => ["Counts"], :type => :single)

  #  tsv.through do |k, values|
  #    tf, tg = values.values_at("TF Associated Gene Name", "TG Associated Gene Name")
  #    pair = [tf, tg] * ":"
  #    triplet = k.split(":").values_at(0,2,3) * ":"
  #    pmid_counts[pair] ||= 0
  #    pmid_counts[pair] += 1
  #    sentence_counts[triplet] ||= 0
  #    sentence_counts[triplet] += 1
  #  end

  #  tsv.add_field "PMID counts" do |k,v|
  #    tf, tg = v.values_at("TF Associated Gene Name", "TG Associated Gene Name")
  #    pair = [tf, tg] * ":"
  #    pmid_counts[pair]
  #  end

  #  tsv.add_field "Sentence counts" do |k,v|
  #    triplet = k.split(":").values_at(0,2,3) * ":"
  #    sentence_counts[triplet]
  #  end

  #  tsv
  #end


end
