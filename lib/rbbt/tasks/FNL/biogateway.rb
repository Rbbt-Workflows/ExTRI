module FNL
  dep :sentence_coverage_NER
  task :biogateway => :tsv do
    tsv = step(:sentence_coverage_NER).load

    tsv = tsv.swap_id "Transcription Factor (Associated Gene Name)", "UniProt/SwissProt Accession"
    tsv.rename_field "UniProt/SwissProt Accession", "Transcription Factor (UniProt/SwissProt Accession)"

    tsv = tsv.swap_id "Target Gene (Associated Gene Name)", "Entrez Gene ID"
    tsv.rename_field "Entrez Gene ID", "Target Gene (Entrez Gene ID)"

    tsv
  end
end
