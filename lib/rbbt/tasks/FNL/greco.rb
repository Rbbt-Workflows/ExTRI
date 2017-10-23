module FNL

  dep :validation_dataset
  dep :FNL_clean
  task :greco_format => :text do 
    tsv = step(:validation_dataset).load.select("Valid" => "Valid")
    fixed = step(:FNL_clean).load

    name2ens = Organism.identifiers("Hsa/feb2014").index :persist => true
    hashes = []
    tsv.through do |pair,values|
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


    hashes.to_json
  end
end
