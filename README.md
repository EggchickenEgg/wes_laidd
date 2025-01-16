# GC% threshold analysis pipeline 1.17
![image](https://github.com/user-attachments/assets/1fa2e34f-32f5-48e9-87ed-ed73c659788c)
import pandas as pd
import subprocess

# find_GC_transition_site

# 설정
input_file = "processed_gc_content.csv"  # GC% 데이터 파일

# 데이터 로드
gc_data = pd.read_csv(input_file)

for i in range(5,16):
    threshold = i  # ΔGC% 임계값

    # ΔGC% 계산
    gc_data["ΔGC%"] = gc_data["GC_Content"].diff().abs()
    gc_data = gc_data[gc_data["ΔGC%"] >= threshold]
    gc_data.to_csv(f"processed_gc_content_threshold_{threshold}", sep = "\t", header = False, index =False)

    print(f'threshold_{threshold}_저장_완료')

### GC_intersect_with_mutation

import os

compare_list = ["data_mutations_leukemia.bed",'data_mutations_melanoma.bed']

results = []

for list in compare_list:
    for i in range(5,16):
                                                
        # 파일 로드
        intersect_file = f"inter_thres_{i}_{list}.bed"
        gc_transition_file = f"processed_gc_content_threshold_{i}"
        compare_file = list

        # if os.path.exists(intersect_file):
        #     print(f"{intersect_file} 파일이 이미 존재합니다. 과정을 건너뜁니다.")
        # else:
        #     print('intersected_file_존재하지 않습니다.')
        #     print('intersected_file을 생성합니다.')
        #     # intersect

        command = f"bedtools intersect -a {compare_file} -b {gc_transition_file} -wa -wb > {intersect_file}"
        result = subprocess.run(command, shell=True)

        if result.returncode == 0:
            print("명령어 실행 성공!")
        else:
            print("명령어 실행 실패:", result.stderr)

        # BED 파일 로드
        gc_data = pd.read_csv(gc_transition_file, sep="\t", header=None)
        compare_data = pd.read_csv(compare_file, sep="\t", header=0)  # 헤더 포함
        intersect_data = pd.read_csv(intersect_file, sep="\t", header=None)

        # 교집합 비율 계산
        total_variants = len(compare_data)
        intersect_variants = len(intersect_data)

        intersection_ratio = (intersect_variants / total_variants) * 100

        # 결과 출력
        print(f"threshold{i}")
        print(f"총 변이 수: {total_variants}")
        print(f"GC% 전이 구역과 일치하는 변이 수: {intersect_variants}")
        print(f"일치 비율: {intersection_ratio:.2f}%")

        #결과 저장
        results.append({
        "Threshold": i,
        "Total Variants": total_variants,
        "Intersect Variants": intersect_variants,
        "Intersection Ratio (%)": f"{intersection_ratio:.2f}",
        "compared_data_type:":list})
        

        # 교집합된 유전자 리스트 저장
        intersected_genes = intersect_data[3].unique()
        pd.DataFrame(intersected_genes, columns=[3]).to_csv(f"{compare_file}_with_{gc_transition_file}.csv", index=False)
        print(f"교집합된 유전자 리스트 저장 완료: intersected_genes_threshold_{i}.csv")


results_df = pd.DataFrame(results)

# 결과 출력
print(results_df)

# CSV 파일로 저장
results_df.to_csv(f"intersection_results_{compare_list[0]}&{compare_list[2]}.csv", index=False)


