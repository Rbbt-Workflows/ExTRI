module FNL


  dep :pairs
  task :all_pmids => :array do
    tsv = step(:pairs).load
    pmid_fields = tsv.fields.select{|f| f=~ /PMID/}
    tsv.slice(pmid_fields).values.flatten.compact.uniq.collect{|v| v.split(/[;|]/)}.flatten.uniq.select{|v| v =~ /^\d{5,9}$/ }.collect{|pmid| pmid.sub(/^0+/,'')}
  end

  dep :pairs
  input :tf_pair, :select, "Use TF or TF-TG pairs", "TF", :select_options => ["TF", "TF-TG"]
  input :databases, :array, "Databases to compare", ["FNL", "TRRUST", "HTRI", "TFacts", "Intact"]
  input :confidence_db, :boolean, "Filter DB entries for high confidence", false
  input :confidence_FNL, :boolean, "Filter FNL entries for high confidence", false
  input :remove_autoregulation, :boolean, "Filter out FNL entries for auto-regulation", false
  extension :png
  task :venn => :binary do |tf_pair,databases,confidence_db,confidence_FNL,remove_autoregulation|
    tsv = step(:pairs).load

    tsv = tsv.select("Auto-regulation"){|v| v.empty?} if remove_autoregulation

    if confidence_db
      good_confidence_fields = tsv.fields.select{|f| f =~ /conf/i}.select{|f| databases.select{|db| f.include? db}.any?}
      good_confidence_fields.each do |f|
        name = f.match(/\[(.*)\]/)[1]
        next if name == "FNL"
        tsv.select(f){|v| v =~ /Low/i}.keys.each do |k|
          tsv[k]["[#{ name }] present"] = ""
        end
      end
    end

    if confidence_FNL
      tsv.select("[FNL] Confidence" => "Low").keys.each do |k|
        tsv[k]["[FNL] present"] = ""
      end
    end

    good_present_fields = tsv.fields.select{|f| f =~ /present/}.select{|f| databases.select{|db| f.include? db}.any?}

    present = nil
    if tf_pair == "TF"
      tsv.with_monitor self.progress_bar("Reordering by TF") do
        present = TSV.setup({}, :key_field => tsv.key_field, :fields => good_present_fields, :type => :list)
        tsv.slice(good_present_fields).through do |k,values|
          k = k.split(":").first
          current = present[k] || values
          new = current.zip(values).collect{|c,n| [c,n].reject{|v| v.nil? or v.empty?}.first }
          present[k] =  new
        end
        Log.tsv present
      end
    else
      present = tsv.slice(good_present_fields)
    end

    Open.write(file('present'), present.to_s)

    log :producing
    require 'rbbt/util/R'
    present.R <<-EOF
    data = apply(data,2,function(x){ x != ""})
    fields = #{R.ruby2R good_present_fields}
    categories = #{R.ruby2R good_present_fields.collect{|f| num = present.select(f){|l| ! (l.nil? || l.empty?) }.length; [f.match(/\[(.*)\]/)[1], num] * ": "}}
    png('#{self.path}')
    rbbt.plot.venn(data, fields, category=categories, margin=0.2)
    dev.off()
    EOF
    nil
  end

  #dep :FNL_confidence
  #input :confidence, :select, "Confidence criteria", "Predicted", :select_options => ["Predicted", "Threshold"]
  #task :top_pairs => :tsv do |confidence|
  #  tsv = step(:FNL_confidence).load

  #  confidence_field = confidence == "Predicted" ? "Prediction confidence" : tsv.fields.select{|f| f =~ /Threshold/}.first
  #  all = TSV.setup({}, :key_field => "TF:TF", :fields => ["FNL all sentences"], :type => :single, :cast => :to_i)
  #  high = TSV.setup({}, :key_field => "TF:TF", :fields => ["FNL High conf. sentences"], :type => :single, :cast => :to_i)
  #  low = TSV.setup({}, :key_field => "TF:TF", :fields => ["FNL Low conf. sentences"], :type => :single, :cast => :to_i)
  #  tsv.through do |key, values|
  #    pair = key.split(":").values_at -2, -1
  #    pair = pair * ":"
  #    conf = values[confidence_field]
  #    all[pair] ||= 0
  #    all[pair] += 1
  #    if conf == "High"
  #      high[pair] ||= 0
  #      high[pair] += 1
  #    else
  #      low[pair] ||= 0
  #      low[pair] += 1
  #    end
  #  end

  #  all.to_list.attach(high.to_list).attach(low.to_list)
  #end
  
  dep :pairs
  dep :triplet_confidence
  input :db, :select, "Database to consider for PMID counts", "FNL", :select_options => FNL::DATABASES
  input :type, :select, "Aggregate by TF, TG, or pair", "TF", :select_options => %w(TF TG TF:TG)
  input :by_triplet, :boolean, "Consider confidence by triplet in FNL", true
  input :remove_autoregulation, :boolean, "Filter out FNL entries for auto-regulation", false
  task :articles => :tsv do |db,type, by_triplet,remove_autoregulation|
    tsv = step(:pairs).load

    tsv = tsv.select("Auto-regulation"){|v| v.empty?} if remove_autoregulation

    triplet_confidence = step(:triplet_confidence).load

    confidence_field = tsv.fields.select{|f| f.include?(db) and f.include?('onfidence')}.first
    pmid_field = tsv.fields.select{|f| f.include?(db) and f.include?('PMID')}.first
    raise ParameterException, "NO PMID field for db: " << db if pmid_field.nil?

    tsv = tsv.slice([pmid_field, confidence_field].compact)

    counts = TSV.setup({}, :key_field => type, :fields => %w(All High Low), :type => :double)
    tsv.through do |k,values|
      tf, tg = k.split(":")
      pmids, conf = values

      
      key = case type
            when "TF"
              tf
            when "TG"
              tg
            else
              k
            end

      next if pmids.nil? or pmids.empty?
      pmids = pmids.split(";")

      counts[key] ||= [[],[],[]]

      counts[key][0].concat pmids

      if db == "FNL" and by_triplet
        pmids_h = []
        pmids_l = []
        pmids.each do |p|
          t = [tf,tg,p] * ":"
          case triplet_confidence[t]
          when "High"
            pmids_h << p
          when "Low"
            pmids_l << p
          else 
            raise "TOP"
          end
        end

        counts[key][1].concat pmids_h
        counts[key][2].concat pmids_l
      else
        if conf == "High"
          counts[key][1].concat pmids
        end

        if conf == "Low"
          counts[key][2].concat pmids
        end
      end

      counts[key].each do |l| l.uniq! end
    end

    counts
  end

  dep :articles
  task :top => :tsv do
    tsv = step(:articles).load
    type = step(:articles).inputs[:type]
    new = tsv.annotate({})
    new.type = :list
    tsv.through do |e,counts|
      new[e] = counts.collect{|l| l.length}
    end

    new.add_field "TF" do |k,v|
      k.split(":").first
    end if type == "TF:TG"
     

    new.add_field "TG" do |k,v|
      k.split(":").last
    end if type == "TF:TG"

    new
  end

  export :venn, :top

  dep :FNL_clean
  task :gene_sentence_counts => :tsv do
    counts = TSV.setup({}, :key_field => "Associated Gene Name", :fields => ["Sentences as TF", "Sentences as TG"], :type => :list, :namespace => FNL.organism)
    TSV.traverse step(:FNL_clean), :bar => self.progress_bar("Counting genes") do |id,values|
      tf, tg, *rest = values
      counts[tf] ||= [0,0]
      counts[tg] ||= [0,0]
      counts[tf] = [counts[tf][0] + 1, counts[tf][1]]
      counts[tg] = [counts[tg][0], counts[tg][1] + 1]
    end

    tsv = step(:FNL_clean).load
    tsv.monitor = true

    word_counts = {}
    tsv.through do |k,values|
      tf, tg, *rest = values
      word_counts[tf] ||= {}
      word_counts[tg] ||= {}
      sentence = values["Sentence"].downcase
      sentence.split(/[^\w]/).each do |w|
        word_counts[tf][w] ||= 0
        word_counts[tf][w] += 1
        word_counts[tg][w] ||= 0
        word_counts[tg][w] += 1
      end
    end

    counts.with_monitor do 
      counts.add_field "Top words" do |gene,values|
        words = word_counts[gene]
        words.sort_by{|w,c| c}.reverse[0..10].collect{|w,c| "#{w} (#{c})" } * ", "
      end
    end

    counts
  end

  dep :FNL_clean
  input :gene, :string, "Gene name"
  task :words => :tsv do |gene|
    tsv = step(:FNL_clean).load
    counts = {} 
    tsv.with_monitor do
      tsv.through do |k,v|
        tf, tg, *rest = v
        next unless tf == gene or tg == gene
        sentence = v["Sentence"].downcase
        sentence.split(/[^\w]/).each do |w|
          counts[w] ||= 0
          counts[w] +=1
        end
      end
    end
    TSV.setup(counts, :key_field => "Word", :fields => ["Counts"], :type => :single, :cast => :to_i)
  end


end
