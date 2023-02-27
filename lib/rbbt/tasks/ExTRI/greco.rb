module ExTRI

  dep :validation_dataset
  dep :ExTRI_clean
  task :greco_format_validation => :text do 
    tsv = step(:validation_dataset).load.select("Valid" => "Valid")
    fixed = step(:ExTRI_clean).load

    name2ens = Organism.identifiers(ExTRI.organism).index :persist => true
    hashes = []
    tsv.through do |pair,values|
      next unless fixed.include? pair
      sentence = fixed[pair]["Sentence"]
      pmid, n, tf, tg = pair.split(":")

      ens_tf = name2ens[tf]
      ens_tg = name2ens[tg]

      #url_tf = "http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=#{ens_tf}"
      #url_tg = "http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=#{ens_tg}"
      #
      url_tf = "http://identifiers.org/ensembl:#{ens_tf}"
      url_tg = "http://identifiers.org/ensembl:#{ens_tg}"
      
      
      hashes << {
        "ext_id": pmid, #pmid
        "provider": "ExTRI", 
        "anns": [
          {
            "exact": sentence,
            "section": "abstract",
            "tags": [
              {"name": tf,
               "uri": url_tf },
               {"name": tg,
                "uri": url_tg }
            ]
          }
        ]
      }
    end


    hashes.to_json
  end

  dep :ExTRI_confidence
  task :europePMC_format => :text do 
    fixed = step(:ExTRI_confidence).load

    hc = fixed.select("Prediction confidence" => "High")

    name2ens = Organism.identifiers(ExTRI.organism).index :persist => true
    hashes = []
    hc.through do |pair,values|
      next unless fixed.include? pair
      sentence = fixed[pair]["Sentence"]
      pmid, n, tf, tg = pair.split(":")

      ens_tf = name2ens[tf]
      ens_tg = name2ens[tg]

      #url_tf = "http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=#{ens_tf}"
      #url_tg = "http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=#{ens_tg}"
      #
      url_tf = "http://identifiers.org/ensembl:#{ens_tf}"
      url_tg = "http://identifiers.org/ensembl:#{ens_tg}"
      
      
      hashes << {
        "src": "MED",
        "id": pmid,
        "provider": "FNL", 
        "anns": [
          {
            "exact": sentence,
            "section": "abstract",
            "tags": [
              {"name": tf,
               "uri": url_tf },
               {"name": tg,
                "uri": url_tg }
            ]
          }
        ]
      }
    end

    Misc.ordered_divide(hashes, 10_000).each_with_index do |chunk,i|
      Open.write(file("ExTRI_HC_#{Time.now.to_date}_chunk_" + i.to_s), chunk.collect{|h| h.to_json } * "\n")
    end


    hashes.collect{|h| h.to_json} * "\n"
  end
end
