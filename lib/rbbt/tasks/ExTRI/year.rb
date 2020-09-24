module ExTRI

  dep :all_pmids
  task :pmid_years => :tsv do
    all_pmids = step(:all_pmids).load
    require 'rbbt/sources/pubmed'

    log :loading_PMIDs
    articles = PubMed.get_article(all_pmids)

    tsv = TSV.setup({}, :key_field => "PMID", :fields => ["Year"], :type => :single)
    log :processing_years
    articles.each do |pmid, article|
      year = article.year
      tsv[pmid] = year
    end

    tsv
  end

  dep :all_pmids
  task :pmid_journal => :tsv do
    all_pmids = step(:all_pmids).load
    require 'rbbt/sources/pubmed'
    articles = PubMed.get_article(all_pmids)

    tsv = TSV.setup({}, :key_field => "PMID", :fields => ["Journal"], :type => :single)
    articles.each do |pmid, article|
      journal = article.journal
      tsv[pmid] = journal
    end

    tsv
  end


  dep :pairs
  dep :pmid_years
  dep :pmid_journal
  input :high_confidence, :boolean, "Consider only High-confidence", false
  task :pmid_details => :tsv do |high_confidence|
    tsv = step(:pairs).load
    pmid_years = step(:pmid_years).load
    pmid_journal = step(:pmid_journal).load
    all_fields = tsv.fields

    databases = ExTRI::DATABASES
    pmid_database = {}
    tsv.with_monitor do
      tsv.through do |pair, values|
        databases.each do |database|
          pmid_field = "[#{database}] PMID"
          next unless all_fields.include? pmid_field
          confidence_field = "[#{database}] Confidence"
          confidence = all_fields.include?(confidence_field)? values[confidence_field] : ""
          next if confidence.include?("Low") && high_confidence
          pmids = values[pmid_field]
          pmids.collect{|p| p.split(";")}.flatten.uniq.each do |pmid|
            pmid = pmid.split(":").last
            pmid_database[pmid] ||= []
            pmid_database[pmid] << [pair, database]
          end
        end
      end
    end

    dumper = TSV::Dumper.new :key_field => "PMID", :fields => ["Year", "Journal"] + databases, :type => :double
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
      year = pmid_years[pmid] 
      journal = pmid_journal[pmid] 
      next if year.nil?
      [pmid, [year, journal] + values]
    end

    dumper
  end

  dep :pmid_details
  task :database_years => :tsv do
    pmid_details = step(:pmid_details).load

    years = pmid_details.column("Year").values.flatten.compact.uniq.sort

    database_years = TSV.setup({}, :key_field => "Database", :fields => years, :type => :double)
    databases = pmid_details.fields[2..-1]
    TSV.traverse pmid_details do |pmid, values|
      year = values["Year"].first
      next if year.nil?
      res = []
      year_counts = {}
      databases.zip(values[2..-1]).each do |database, values|
        if values.length > 0
          database_years[database] ||= [0] * years.length
          database_years[database][years.index(year)] += 1
        end
      end
    end
    database_years
  end

  dep :database_years
  extension :png
  task :database_year_plot => :binary do
    log :producing
    require 'rbbt/util/R'

    step(:database_years).load.R <<-EOF
rbbt.require('reshape')
rbbt.require('plyr')
rbbt.require('dplyr')
rbbt.require('ggplot2')

m <- rbbt.tsv.melt(data, "Year", "Count")
m <- subset(m, Count > 0)
remove.GEREDB <- subset(m, Database == "GEREDB" & Year==2011 )
remove.ExTRI <- subset(m, Database == "ExTRI" & Year==2014 )
m = anti_join(m, remove.GEREDB)
m = anti_join(m, remove.ExTRI)

years = as.integer(levels(m$Year))
breaks = pretty(years, length(years) / 2)
g <- ggplot(m, aes(x=Year, y=Count)) + geom_line(aes(group=Database, color=Database, linetype=Database)) + scale_x_discrete(breaks=breaks) + 
   theme_classic() + 
   theme(
        axis.title.y  = element_text(angle=90, vjust=0.5, size=20), 
        axis.text.y  = element_text(angle=90, vjust=0.5, size=20), 
        axis.title.x  = element_text(size=20), 
        axis.text.x  = element_text(angle=90, vjust=0.5, size=20), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=16))

rbbt.png_plot('#{self.tmp_path}', 'plot(g)', height=600, width=800)

data = NULL
    EOF

    nil
  end

  dep :pmid_details
  task :TF_years => :tsv do
    pmid_details = step(:pmid_details).load

    tf_years = TSV::Dumper.new(:key_field => "TF", :fields => ["Years"], :type => :flat)
    tf_years.init
    TSV.traverse pmid_details, :into => tf_years do |pmid, values|
      year, journal, *rest = values
      res = rest.flatten.collect{|p| p.split(":").first}.uniq.collect{|tf| [tf, year]}
      res.extend MultipleResult
      res
    end
    
    TSV.collapse_stream(tf_years.stream, :uniq => false)
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

  dep :TF_years
  input :genes, :array, "Genes to consider (empty to consider all)"
  input :top, :integer, "Show only top genes", 100
  input :equal_height, :boolean, "Show bar with equal height", false
  input :remove_autoregulation, :boolean, "Filter out ExTRI entries for auto-regulation", false
  input :remove_non_TFClass, :boolean, "Filter out ExTRI entries for non TFClass TF", false
  extension :svg
  task :life_cycle => :text do |genes,top,equal_height,remove_autoregulation,remove_non_TFClass|
    tsv = step(:TF_years).load
    tsv = tsv.select(genes) if genes and genes.any?

    tsv = tsv.select("Auto-regulation"){|v| v.empty?} if remove_autoregulation
    tsv = tsv.select("TFClass"){|v| ! v.empty?} if remove_non_TFClass

    counts = {}

    tsv.through do |tf, articles|
      next unless articles.length > 10
      counts[tf.gsub('-', '_')] = Misc.counts(articles.collect{|pmid| "y" << pmid.gsub("-", '_')})
    end

    require 'rbbt/util/R'
    log :producing
    R.run  <<-EOF, [:svg]
      library(reshape)
      library(ggplot2)

      data = #{R.ruby2R(counts).gsub("),", "),\n")}

      m = melt(data)
      names(m) <- c("value", "Year", "Gene")

      m = subset(m, m$Year!='y')

      gene.counts = aggregate(value ~ Gene, m, sum)
      rownames(gene.counts) <- gene.counts$Gene

      m$Gene = factor(m$Gene, levels = unique(m$Gene[sort(gene.counts[m$Gene,'value'], index.return=T, decreasing=T)$ix]))

      m = subset(m, m$Gene %in% rownames(gene.counts)[1:#{top}])

      g <- ggplot(m, aes(x=Gene, fill=Year, weight=value)) 

      #{equal_height ?  "g <- g + geom_bar(position='fill')": "g <- g + geom_bar()" }
      
      g <- g + theme_gdocs() + theme(axis.text.x = element_text(angle = 60, hjust=1)) 

      rbbt.SVG.save(file='#{self.path}', g, width=15, height=5)
    EOF

    nil
  end
end
