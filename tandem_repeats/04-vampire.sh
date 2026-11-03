max_jobs=20
log_dir="run_logs"
mkdir -p "$log_dir"

for dir in ./*/; do
    [ -d "$dir" ] || continue
    name=$(basename "$dir")
    out_prefix="${name}_TR_anno"
    log_file="${log_dir}/${name}.log"
    cd ${dir}
    # 后台执行，输出日志，任务失败直接exit 1标记
    vampire anno \
        --no-auto \
        --no-denovo \
        "all.sample.fa" \
        "$out_prefix" \
        -m "${name}.fa" \
        >"$log_file" 2>&1 || exit 1 &
    cd ..
    if [ $(jobs -p | wc -l) -ge "$max_jobs" ]; then
        wait -n
    fi
done

wait
echo "所有任务执行完毕，日志存放：$log_dir"
