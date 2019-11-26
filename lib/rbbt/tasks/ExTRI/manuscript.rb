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
  task :table_pairs_content => :tsv do
    data = step(:pairs).load
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

    res
  end
end
