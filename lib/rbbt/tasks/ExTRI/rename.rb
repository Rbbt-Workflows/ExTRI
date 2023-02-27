module ExTRI
  def self.translations
    @@translations ||= Rbbt.root.data.update_gene_names_v2.tsv
  end

  def self.update_symbols(value)
    if Array === value
      value.collect{|v| update_symbols(v) }
    else
      updated = translations[value]
      Log.debug "Updated gene symbol #{value} => #{updated}" if updated
      updated || value
    end
  end

  def self.update_key_symbols(key)
    update_symbols(key.split(":")) * ":"
  end

  def self.update_tsv_symbols(old)
    old = TSV.open(old) unless TSV === old
    new = old.annotate({})
    TSV.traverse old, :into => new do |k,v|
      k = ExTRI.update_key_symbols(k)
      v = ExTRI.update_symbols(v)

      [k, v]
    end
  end

  dep :pairs_final
  task :gene_names => :array do
    step(:pairs_final).load.keys.collect{|k| k.split(":") }.flatten.uniq.sort
  end


  dep :gene_names
  task :bad_names => :tsv do
    good_names = Organism.identifiers(ExTRI.updated_organism).tsv(:key_field => "Associated Gene Name", :fields => []).keys

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
    index = Organism.identifiers(ExTRI.updated_organism).index :target => "Associated Gene Name", :persist => true, :order => true

    dumper = TSV::Dumper.new(:key_field => "original (Associated Gene Name)", :fields => ["TF", "good (Associated Gene Name)"], :type => :list)
    dumper.init
    TSV.traverse step(:bad_names), :type => :single, :into => dumper do |gene,tf|
      good = index[gene]
      [gene, [tf, good]]
    end
  end

  input :old_tsv, :tsv, "Old tsv"
  task :update_tsv_symbols => :tsv do |old|
    ExTRI.update_tsv_symbols(old)
  end

end
