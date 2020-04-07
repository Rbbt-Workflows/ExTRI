require 'rbbt/util/R'
module ExTRI

  def self.count_nodes(info, count = {}, good = nil)
    name = info["name"] 

    if info["children"].nil?
      return [[], 0] if info["id"].split(".").length < 5
      parts = name.scan(/[a-z0-9\-\?]+/i)
      parts = parts & good if good
      return [parts, count] 
    end

    total = []
    info["children"].each do |child|
      list, _c = count_nodes(child, count, good)
      total += list
    end

    count[name] = total
    [total, count]
  end

  dep :ExTRI_confidence
  dep :top, :type => "TF"
  input :high_confidence, :boolean, "Restrict to high confidence predictions", false
  task :sunburst_percents => :yaml do |high_confidence|
    tsv = step(:ExTRI_confidence).load
    top = step(:top).load

    tsv = tsv.select("Prediction confidence" => 'High') if high_confidence
    tfs = tsv.column("Transcription Factor (Associated Gene Name)").values.uniq

    tree = JSON.parse(TFClass.hierarchy_json.produce.read)
    baseline = ExTRI.count_nodes(tree).last
    good = ExTRI.count_nodes(tree, {}, tfs).last

    good_counts = TSV.setup({}, "Catergory~All,Found#:type=:list#:cast=:to_f")
    good.each do |cat,genes|
      all = baseline[cat]
      good_counts[cat] = [all.length, genes.length]
    end

    Open.write(file('genes_per_cat.yaml'), good.to_yaml)
    Open.write(file('counts_per_cat.tsv'), good_counts.to_s)

    percent = {}
    baseline.each do |k,v|
      genes = good[k]
      p = genes.any? ? genes.length.to_f / v.length : 0
      top_p = top.select(genes).sort_by{|g,c| c.first}.reverse.collect{|g,c| g}[0..2]

      case top.select(genes).length
      when 0
      when 1, 2, 3
        k = k + " (e.g. " + (top_p * ", ") + ")" if top_p.any?
      else
        k = k + " (e.g. " + (top_p * ", ") + ", ...)" if top_p.any?
      end

      percent[k] = p
    end

    percent
  end

  dep :pairs
  dep :top, :type => "TF", :db => "All KB"
  input :remove_ExTRI, :boolean, "Do not include ExTRI in the counts when looking at All KB", false
  input :high_confidence, :boolean, "Restrict ExTRI to high confidence predictions", false
  task :sunburst_percents_kb => :yaml do |remove_ExTRI,high_confidence|
    tsv = step(:pairs).load
    top = step(:top).load

    if remove_ExTRI
      presence_fields = tsv.fields.select{|f| f.include?('present') && ! f.include?('ExTRI')} 
    else
      presence_fields = tsv.fields.select{|f| f.include?('present')} 
    end
    
    tfs = presence_fields.inject([]) do |acc,f|
      if f.include?("ExTRI") and high_confidence
        acc += tsv.select("[ExTRI] Confidence" => 'High').select(f){|v| ! (v.empty? || v.nil?)  }.keys.collect{|k| k.split(":").first}
      else
        acc += tsv.select(f){|v| ! (v.empty? || v.nil?)  }.keys.collect{|k| k.split(":").first}
      end
    end.uniq
    
    tree = JSON.parse(TFClass.hierarchy_json.produce.read)
    baseline = ExTRI.count_nodes(tree).last
    good = ExTRI.count_nodes(tree, {}, tfs).last

    good_counts = TSV.setup({}, "Catergory~All,Found#:type=:list#:cast=:to_f")
    good.each do |cat,genes|
      all = baseline[cat]
      good_counts[cat] = [all.length, genes.length]
    end

    Open.write(file('genes_per_cat.yaml'), good.to_yaml)
    Open.write(file('counts_per_cat.tsv'), good_counts.to_s)
    percent = {}
    baseline.each do |k,v|
      genes = good[k]
      p = genes.any? ? genes.length.to_f / v.length : 0
      top_p = top.select(genes).sort_by{|g,c| c.first}.reverse.collect{|g,c| g}[0..2]
      case top.select(genes).length
      when 0
      when 1, 2, 3
        k = k + " (e.g. " + (top_p * ", ") + ")" if top_p.any?
      else
        k = k + " (e.g. " + (top_p * ", ") + ", ...)" if top_p.any?
      end
      percent[k] = p
    end

    percent
  end

  dep :sunburst_percents_kb do |jobname,options|
    if options[:source].to_s == "ExTRI"
      {:task => :sunburst_percents, :jobname => jobname, :inputs => options}
    else
      {:task => :sunburst_percents_kb, :jobname => jobname, :inputs => options}
    end
  end
  input :source, :select, "Source of TRI", :ExTRI, :select_options => %w(ExTRI KB)
  task :sunburst => :text do 
    percent = dependencies.first.load

    labels = percent.keys
    percents = percent.values_at *labels

    tree_text = TFClass.hierarchy_json.produce.read

    labels.each do |label|
      orig = label.split(" (e.g.").first
      tree_text.gsub!(/"#{orig}"/, '"' + label + '"')
    end

    new_tree = {}

    Open.write(file('tree.json'), tree_text)
    
    widget = file('widget.html')
    Open.mkdir self.files_dir
    TmpFile.with_file(tree_text) do |tree_file|
      TmpFile.with_file(labels.to_json) do |labels|
        TmpFile.with_file(percents.to_json) do |percents|
          R.run <<-EOF
rbbt.require('timelyportfolio/sunburstR')
rbbt.require("rjson")
rbbt.require("colorspace")

labels = fromJSON(file='#{labels}')
percents = fromJSON(file='#{percents}')
#basic_colors = rep(rev(RColorBrewer::brewer.pal(11, "BrBG")), length(labels)/11)
basic_colors = rep("#000000", length(labels))

colors = c()
white = RGB(255,255,255)
for (i in seq(1,length(labels))){
  percent = percents[i]
  basic_color = col2rgb(basic_colors[i])
  color = RGB(basic_color["red",1], basic_color["green",1], basic_color["blue",1])
  new_color_rgb = mixcolor(percent, white, color)@coords[1,]
  new_color = rgb(new_color_rgb[1], new_color_rgb[2], new_color_rgb[3],max=255)

  colors = c(colors, new_color)
}

rbbt.require('readr')
data=read_file('#{tree_file}')

#sunburst(data, colors=list(range=colors,domain=labels))

options = list(legendOrder = NULL, colors = list(range=colors,domain=labels),
valueField = "size", percent = TRUE, count = FALSE, explanation = NULL,
breadcrumb = list(), legend = list(), sortFunction = NULL, sumNodes = TRUE,
withD3 = FALSE, width = NULL, height = NULL, elementId = NULL, sizingPolicy =
NULL, csvdata = NULL, jsondata = NULL)

data.fixed = jsonlite::toJSON(jsonlite::fromJSON(data), auto_unbox = TRUE, dataframe = "rows") 

x = list(data=data.fixed, options=options)
width = NULL
height = NULL
elementId = NULL
dep <- d3r::d3_dep_v4()
sizingPolicy <- htmlwidgets::sizingPolicy(browser.fill = TRUE)

w <- htmlwidgets::createWidget(name = "sunburst", x, width = width,
height = height, package = "sunburstR", elementId = elementId,
sizingPolicy = sizingPolicy, dependencies = dep)

htmlwidgets::saveWidget(w, file='#{widget}', selfcontained=FALSE)

EOF
        end
      end
    end
    "DONE"
  end
end
