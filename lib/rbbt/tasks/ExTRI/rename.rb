module ExTRI
  def self.translations
    @@translations ||= Rbbt.root.data.update_gene_names.tsv
  end

  def self.fix_name(value)
    if Array === value
      value.collect{|v| fix_name(v) }
    else
      translations[value] || value
    end
  end

  def self.fix_tsv_names(old)
    old = TSV.open(old) unless TSV === old
    new = old.annotate({})
    TSV.traverse old, :into => new do |k,v|
      k = ExTRI.fix_name(k.split(":")) * ":"

      [k, ExTRI.fix_name(v)]
    end
  end

  dep :pairs_final
  task :gene_names => :array do
    step(:pairs_final).load.keys.collect{|k| k.split(":") }.flatten.uniq.sort
  end


  dep :gene_names
  task :bad_names => :tsv do
    good_names = Organism.identifiers("Hsa/feb2021").tsv(:key_field => "Associated Gene Name", :fields => []).keys

    tfs = step(:pairs_final).load.keys.collect{|k| k.split(":").first }
    bad_names = step(:gene_names).load - good_names

    tsv = TSV.setup(bad_names, :key_field => "Associated Gene Name", :fields => [])
    tsv.add_field "TF" do |k,v|
      tfs.include? k
    end
    tsv
  end


  dep :bad_names
  task :good_names => :tsv do
    index = Organism.identifiers("Hsa/feb2021").index :target => "Associated Gene Name", :persist => true, :order => true

    dumper = TSV::Dumper.new(:key_field => "original (Associated Gene Name)", :fields => ["TF", "good (Associated Gene Name)"], :type => :single)
    dumper.init
    TSV.traverse step(:bad_names), :type => :single, :into => dumper do |gene,tf|
      good = index[gene]
      [gene, [tf, good]]
    end
  end

  input :old_tsv, :tsv, "Old tsv"
  task :fix_tsv_names => :tsv do |old|
    ExTRI.fix_tsv_names(old)
  end

end
