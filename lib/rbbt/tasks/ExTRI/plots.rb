require 'rbbt/util/R'
module ExTRI

  def self.count_nodes(info, count = {}, good = nil)
    name = info["name"] 
    parts = name.split(/[\s,]+/)
    total = (good.nil? || (good & parts).any?) ? 1 : 0
    return [total, count] unless info["children"]
    info["children"].each do |child|
      sub_total, _c = count_nodes(child, count, good)
      total += sub_total
    end
    count[name] = total
    [total, count]
  end

  dep :ExTRI_confidence
  input :high_confidence, :boolean, "Restrict to high confidence predictions", false
  task :sunburst_percents => :yaml do |high_confidence|
    tsv = step(:ExTRI_confidence).load
    tsv = tsv.select("Prediction confidence" => 'High') if high_confidence
    tfs = tsv.column("Transcription Factor (Associated Gene Name)").values.uniq
    tree = JSON.parse(TFClass.hierarchy_json.produce.read)
    baseline = ExTRI.count_nodes(tree).last
    good = ExTRI.count_nodes(tree, {}, tfs).last

    percent = {}
    baseline.each do |k,v|
      p = good[k].to_f / v
      percent[k] = p
    end

    percent
  end

  dep :sunburst_percents
  task :sunburst => :text do 
    percent = step(:sunburst_percents).load

    labels = percent.keys
    percents = percent.values_at *labels
    
    widget = file('widget.html')
    Open.mkdir self.files_dir
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
data=read_file('#{TFClass.hierarchy_json.produce.find}')

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
    "DONE"
  end
end
