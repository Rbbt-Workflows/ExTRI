module FNL

  dep :pairs
  dep :all_pmids
  input :high_confidence, :boolean, "Consider only High-confidence", false
  task :pmid_details => :tsv do |high_confidence|
    tsv = step(:pairs).load
    all_pmids = step(:all_pmids).load
    all_fields = tsv.fields

    databases = %w(FNL HTRI TRRUST TFacts Intact)
    pmid_database = {}
    tsv.with_monitor do
      tsv.through do |pair, values|
        databases.each do |database|
          pmid_field = "[#{database}] PMID"
          confidence_field = "[#{database}] Confidence"
          confidence = all_fields.include?(confidence_field)? values[confidence_field] : ""
          next if confidence == "Low" and high_confidence
          pmids = values[pmid_field]
          pmids.split(";").each do |pmid|
            pmid = pmid.split(":").last
            pmid_database[pmid] ||= []
            pmid_database[pmid] << [pair, database]
          end
        end
      end
    end

    years = Rbbt.data["pmids2year_all.txt"].tsv

    dumper = TSV::Dumper.new :key_field => "PMID", :fields => ["Year"] + databases, :type => :double
    dumper.init
    TSV.traverse pmid_database, :bar => true, :into => dumper do |pmid, entries|
      db_pairs = {}
      databases.each do |database|
        db_pairs[database] ||= []
      end
      entries.each do |pair, database|
        db_pairs[database] << pair
      end
      values = db_pairs.values_at(*databases).collect{|v| v.uniq }
      year = years[pmid] 
      next if year.nil?
      [pmid, [year] + values]
    end

    dumper
  end

  dep :pmid_details
  task :database_years => :tsv do
    pmid_details = step(:pmid_details).load

    database_years = TSV::Dumper.new(:key_field => "Database", :fields => ["Years"], :type => :single)
    databases = pmid_details.fields[1..-1]
    TSV.traverse pmid_details, :into => database_years do |pmid, values|
      year = values.first
      res = []
      databases.zip(values[1..-1]).each do |database, values|
        if values.length > 0
          res << [database, year]
        end
      end
      res.extend MultipleResult
      res
    end
  end

  dep :pmid_details
  task :TF_years => :tsv do
    pmid_details = step(:pmid_details).load

    tf_years = TSV::Dumper.new(:key_field => "TF", :fields => ["Years"], :type => :single)
    TSV.traverse pmid_details, :into => tf_years do |pmid, values|
      year = values.first
      res = values[1..-1].flatten.collect{|p| p.split(":").first}.uniq.collect{|tf| [tf, year]}
      res.extend MultipleResult
      res
    end
  end

  dep :TF_years
  dep :pairs
  task :TF_earliest => :tsv do
    tsv = step(:TF_years).path.tsv :merge => true, :type => :flat
    min = TSV.setup({}, :key_field => "TF", :fields => ["Earliest year"], :type => :single)
    tsv.each do |tf, years|
      min[tf] = years.collect{|y| y.to_i}.min
    end
    min
  end
end
