module ExTRI

  dep :ExTRI_confidence
  task :table_content => :tsv do
    data = step(:ExTRI_confidence).load
    hc_data = data.select("Prediction confidence" => "High")
    res = TSV.setup({},"ExTRI~All,High conf.#:type=:list")

    res["TRIs"] = [data.keys.collect{|k| k.split(":")[2..3]*":"}.uniq.length, hc_data.keys.collect{|k| k.split(":")[2..3]*":"}.uniq.length]
    res["TFs"] = [data.keys.collect{|k| k.split(":")[2]}.uniq.length, hc_data.keys.collect{|k| k.split(":")[2]}.uniq.length]
    res["TGs"] = [data.keys.collect{|k| k.split(":")[3]}.uniq.length, hc_data.keys.collect{|k| k.split(":")[3]}.uniq.length]
    res["Sentences"] = [data.size, hc_data.size]
    res["Abstacts"] = [data.keys.collect{|k| k.split(":")[0]}.uniq.length, hc_data.keys.collect{|k| k.split(":")[0]}.uniq.length]
    res

  end

  dep :pairs
  input :restrict_TFClass, :boolean, "Consider only TFClass", true
  task :table_pairs_content => :tsv do |restrict_TFClass|
    data = step(:pairs).load
    Log.tsv data
    data = data.select("TFClass" => "TFClass") if restrict_TFClass
    Log.tsv data
    res = TSV.setup({},"Database~TF:TG All,TF:TG HC,TF All,TF HC,TG All,TG HC#:type=:list")
    db_pairs = {}
    hc_db_pairs = {}
    ExTRI::DATABASES.each do |db|
      begin
        db_data = data.select("[#{db}] present" => db)
      rescue
        next
      end
      hc_db_pairs[db] = db_pairs[db] = db_data.keys
      begin
        hc_db_pairs[db] = db_data.select("[#{db}] Confidence" => "High").keys
      rescue
      end
    end

    db_pairs_uniq = {}
    hc_db_pairs_uniq = {}
    db_tf_uniq = {}
    db_tg_uniq = {}
    hc_db_tf_uniq = {}
    hc_db_tg_uniq = {}
    all_dbs = db_pairs.keys
    all_dbs.each do |db|
      rest = db_pairs.values_at(*(all_dbs - [db])).flatten.uniq
      hc_rest = hc_db_pairs.values_at(*(all_dbs - [db])).flatten.uniq
      db_pairs_uniq[db] = db_pairs[db] - rest
      hc_db_pairs_uniq[db] = hc_db_pairs[db] - hc_rest

      rest_tf = db_pairs.values_at(*(all_dbs - [db])).flatten.uniq.collect{|p| p.split(":").first}.uniq
      db_tf_uniq[db] = db_pairs[db].collect{|p| p.split(":").first}.uniq - rest_tf

      hc_rest_tf = hc_db_pairs.values_at(*(all_dbs - [db])).flatten.uniq.collect{|p| p.split(":").first}.uniq
      hc_db_tf_uniq[db] = hc_db_pairs[db].collect{|p| p.split(":").first}.uniq - hc_rest_tf

      rest_tg = db_pairs.values_at(*(all_dbs - [db])).flatten.uniq.collect{|p| p.split(":").first}.uniq
      db_tg_uniq[db] = db_pairs[db].collect{|p| p.split(":").first}.uniq - rest_tg

      hc_rest_tg = hc_db_pairs.values_at(*(all_dbs - [db])).flatten.uniq.collect{|p| p.split(":").last}.uniq
      hc_db_tg_uniq[db] = hc_db_pairs[db].collect{|p| p.split(":").last}.uniq - hc_rest_tg


    end

    db_pairs.keys.each do |db|
      res[db] = [
        db_pairs[db],
        hc_db_pairs[db], 
        db_pairs[db].collect{|p| p.split(":").first}.uniq, 
        hc_db_pairs[db].collect{|p| p.split(":").first}.uniq,
        db_pairs[db].collect{|p| p.split(":").last}.uniq, 
        hc_db_pairs[db].collect{|p| p.split(":").last}.uniq,
      ].collect{|l| l.length}

      res[db + "-U"] = [
        db_pairs_uniq[db],
        hc_db_pairs_uniq[db], 
        db_tf_uniq[db], 
        hc_db_tf_uniq[db], 
        db_tg_uniq[db], 
        hc_db_tg_uniq[db], 
      ].collect{|l| l.length}
    end

      res["Union"] = [
        db_pairs.values.flatten.uniq,
        hc_db_pairs.values.flatten.uniq, 
        db_pairs.values.flatten.uniq.collect{|p| p.split(":").first}.uniq, 
        hc_db_pairs.values.flatten.uniq.collect{|p| p.split(":").first}.uniq,
        db_pairs.values.flatten.uniq.collect{|p| p.split(":").last}.uniq, 
        hc_db_pairs.values.flatten.uniq.collect{|p| p.split(":").last}.uniq,
      ].collect{|l| l.length}
    res
  end

  dep :ExTRI_postprocess
  extension :xlsx
  task :post_processing_sentences => :binary do
    book = RubyXL::Workbook.new
    Dir.glob(step(:ExTRI_postprocess).files_dir + "/*").sort.each_with_index do |file,i|
      name = File.basename file
      log :processing, "Processing #{ name }"
      tsv = TSV.open(file)

      fields, rows = TSV._excel_data(tsv)
      
      fields = %w(ID TF TG Score Sentence)

      if i == 0
        sheet = book.worksheets.first
        sheet.sheet_name = name
      else
        sheet = book.add_worksheet(name)
      end

      fields.each_with_index do |e,i|
        sheet.add_cell(0, i, e)
      end

      rows.each_with_index do |cells,i|
        cells.each_with_index do |e,j|
          sheet.add_cell(i+1, j, e)
        end
      end

    end
    book.write self.tmp_path
    nil
  end

  dep :pairs
  task :adhoc => :integer do
    pairs = step(:pairs).load

    selected = pairs.select("[ExTRI] present" => "ExTRI").select("[ExTRI] Confidence" => "High")
    set_info :ExTRI_HC, selected.length
    set_info :ExTRI_HC_low, selected.select("[ExTRI] PMID"){|ls| ls.collect{|l| l.split(";")}.flatten.uniq.length <= 2}.length
    set_info :ExTRI_HC_high, selected.select("[ExTRI] PMID"){|ls| ls.collect{|l| l.split(";")}.flatten.uniq.length > 2}.length
    
    selected.select("[ExTRI] PMID"){|ls| ls.collect{|l| l.split(";")}.flatten.uniq.length <= 2}.select do |k,values|
      values.to_hash.select do |field,value|
        next false if field.include? "ExTRI"
        field.include?('present') && ! value.empty?
      end.any?
    end.length
  end
end
