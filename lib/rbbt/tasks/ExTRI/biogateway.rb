module ExTRI
  #dep :sentence_coverage_NER
  dep :pairs
  task :biogateway => :tsv do
    #tsv = step(:sentence_coverage_NER).load
    tsv = step(:pairs).load

    log :prepare, "Preparing tsv for biogateway translation"
    tsv = tsv.to_double

    log :fixing, "Fixing complexes and families"
    #tsv.process "Transcription Factor (Associated Gene Name)" do |name|
    tsv.process "TF" do |name|
      name = name.first
      case name
      when "AP1"
        AP1_SYN
      when "NFKB"
        NFKB_SYN
      else
        [name]
      end 
    end

    #tsv.process "Target Gene (Associated Gene Name)" do |name|
    tsv.process "TG" do |name|
      name = name.first
      case name
      when "AP1"
        AP1_SYN
      when "NFKB"
        NFKB_SYN
      else
        [name]
      end 
    end


    log :translating, "Translating TF"
    #tsv = tsv.swap_id "Transcription Factor (Associated Gene Name)", "UniProt/SwissProt Accession"
    tsv = tsv.swap_id "TF", "UniProt/SwissProt Accession"
    tsv.rename_field "UniProt/SwissProt Accession", "Transcription Factor (UniProt/SwissProt Accession)"

    log :translating, "Translating TG"
    #tsv = tsv.swap_id "Target Gene (Associated Gene Name)", "Entrez Gene ID"
    tsv = tsv.swap_id "TG", "Entrez Gene ID"
    tsv.rename_field "Entrez Gene ID", "Target Gene (Entrez Gene ID)"

    tsv
  end

  dep :biogateway
  task :biogateway_ids => :tsv do
    tsv = step(:biogateway).load

    ids =<<-EOF.split("\n")
Ensembl Transcript ID
Associated Gene Name
UniProt/SwissProt Accession
Entrez Gene ID
Ensembl Gene ID
Ensembl Protein ID
RefSeq ID
RefSeq mRNA
RefSeq Protein ID
EOF
     
    organism = ExTRI.organism

    ids.each do |id|
      log :add_id, id

      if id != "UniProt/SwissProt Accession"
        index = case id
                when "Ensembl Transcript ID"
                  i = Organism.identifiers(organism).tsv :key_field => "UniProt/SwissProt Accession", :fields => ["Ensembl Gene ID"], :type => :double, :persist => false
                  trans = Organism.transcripts(organism).tsv(:fields => ["Ensembl Transcript ID"], :key_field => "Ensembl Gene ID", :type => :flat, :merge => true, :persist => true)
                  i = i.attach trans
                  i.slice("Ensembl Transcript ID").to_flat
                else
                  Organism.identifiers(organism).tsv :key_field => "UniProt/SwissProt Accession", :fields => [id], :type => :flat, :persist => true
                end

        tsv.add_field "Transcription Factor (#{ id })" do |k,values|
          names = values["Transcription Factor (UniProt/SwissProt Accession)"]
          nids = index.values_at(*names).compact.flatten
          nids
        end
      end

      if id != "Entrez Gene ID"
        index = case id
                when "Ensembl Transcript ID"
                  i = Organism.identifiers(organism).tsv :key_field => "Entrez Gene ID", :fields => ["Ensembl Gene ID"], :type => :double, :persist => false
                  i = i.attach Organism.transcripts(organism), :fields => ["Ensembl Transcript ID"], :persist_input => true
                  i.slice("Ensembl Transcript ID").to_flat
                else
                  Organism.identifiers(organism).tsv :key_field => "Entrez Gene ID", :fields => [id], :type => :flat, :persist => true
                end

        tsv.add_field "Target Gene (#{ id })" do |k,values|
          names = values["Target Gene (Entrez Gene ID)"]
          nids = index.values_at(*names).compact.flatten
          nids
        end
      end
    end

    tsv
  end
end
