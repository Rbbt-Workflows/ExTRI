require 'rbbt/util/R'
module ExTRI

  def self.tree2json(tree)
    tree.collect do |name,subtree|
      next if name == :count
      info = {:name => name, :children => tree2json(subtree)}
      info[:size] = subtree[:count] if subtree[:count]
      info
    end.compact
  end

  extension :svg
  task :sunburst_text => :text do
    txt =<<-EOF
a-b-c
a-b-c-d
a-b-c-d
a-b-c-d
a-b-c-d
a-b-c-d
a-b-c-d
a-b-c-d
a-b-c-d
a-b-c-d
a-b-c-e
a-c
    EOF

    tree = {}
    txt.split("\n").each do |line|
      current = tree
      parts = line.split("-")
      parts.each do |part,i|
        current = current[part] ||= {}
      end
      current[:count] ||= 0
      current[:count] += 1
    end

    data = {:name => 'root', :children => ExTRI.tree2json(tree)}

    R.run <<-EOF
rbbt.require('timelyportfolio/sunburstR')
data = '#{data.to_json}'
sunburst(data)
    EOF
    nil
  end

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
  task :sunburst_percents => :yaml do
    tfs = step(:ExTRI_confidence).load.column("Transcription Factor (Associated Gene Name)").values.uniq
    tree = JSON.parse(Rbbt.data["TFClass.json"].read)
    baseline = ExTRI.count_nodes(tree).last
    good = ExTRI.count_nodes(tree, {}, tfs).last

    percent = {}
    baseline.each do |k,v|
      percent[k] = good[k].to_f / v
    end

    percent
  end

  dep :sunburst_percents
  task :sunburst => :text do
    percent = step(:sunburst_percents).load

    labels = percent.keys
    percents = percent.values_at *labels
    
    TmpFile.with_file(labels.to_json) do |labels|
      TmpFile.with_file(percents.to_json) do |percents|
        R.run <<-EOF
rbbt.require('timelyportfolio/sunburstR')
rbbt.require("rjson")
rbbt.require("colorspace")

labels = fromJSON(file='#{labels}')
percents = fromJSON(file='#{percents}')
basic_colors = rep(rev(RColorBrewer::brewer.pal(11, "BrBG")), length(labels)/11)

colors = c()
white = RGB(0,0,0)
for (i in seq(1,length(labels))){
  percent = percents[i]
  basic_color = col2rgb(basic_colors[i])
  color = RGB(basic_color["red",1], basic_color["green",1], basic_color["blue",1])
  new_color = mixcolor(percent, color, white)
  colors = c(colors, new_color[1,])
}

rbbt.require('readr')
data=read_file('#{Rbbt.data["TFClass.json"].find}')
sunburst(data, colors=list(range=colors,domain=labels))

EOF
      end
    end
  end
end
