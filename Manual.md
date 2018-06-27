# Package 명령어

## 주의
- `from pkprocess import *` 명령어 실행 후에 아래 명령들 실행 가능
- `[key=default]`: optional arguments
- `sd`: SeismicData object
- 헤더 정보는 모든 trace에서 추출
- 헤더 keyword는 [SU Trace header 키워드](http://www.cwp.mines.edu/sututor/node149.html) 참고

## 파일 입출력
- SU 파일 읽기: `sd = read_su("input.su")`
- SU 파일로 쓰기: `sd.write_su("output.su")`
- SeismicData(패키지 자체에서 사용) 읽기: `sd = read("input.sd")`
- SeismicData 쓰기: `sd.write("output.sd")`

## 헤더 정보
- 키 값 얻기(sugethw): `list = get_key(sd, keyword)`
- 여러 개의 키 값 얻기(sugethw): `list = get_keys(sd, list_of_keywords)`
- 키 값 얻기(중복 없이): `list = get_key_unique(sd, keyword)`
- 중복되지 않는 키 값의 개수(sucountkey): `ncount = get_key_count(sd,keyword)`

- 공통송신원모음 개수(keyword 기준으로 송신원 구분): `nshot = get_nshot(sd[, keyword="fldr"])` 
- 트레이스 샘플 개수: `ns = get_ns(sd)`
- 트레이스 샘플링 간격 (초 단위): `dt = get_dt(sd)`

- 전체 트레이스 개수: `ntrace = get_ntr(sd)`
- 각각의 공통송신원에 포함된 트레이스 개수: `ntrace_per_shotgather = ntr_per_shot(sd[, keyword="fldr"])`

- list에 포함된 키 값에 해당하는 트레이스만 추출(suwind): `result_sd = window(sd, keyword, list_of_key_values)`
- 각 키 값의 범위 출력(surange): `print_range(sd)`
- 자료 처리 로그 출력: `sd.print_log([nmo=False|True])`

- 특정 keyword의 같은 키 값을 가지는 트레이스들끼리 모으기(susplit): `list_of_sd = trace_split(sd, keyword)`

## Plotting

- percent 기준으로 데이터 clip: `clipped_sd = perc_clip(sd[, perc=100])`

- 이미지 보기(suximage): `plot_image(sd[, figsize=[5,10], perc=100, cmap="gray_r", ratio=2, f2=0, d2=1, key=False])`
- Wiggle trace 보기(suwigb): `plot_wiggle(data[, figsize=[5,10], fill=True, perc=100, scale=1])`

- 주파수 스펙트럼: `specfx(sd[, perc=100, cmap="jet"])`
- FK 스펙트럼: `specfk(sd[, dx=from_offset, perc=100, cmap="jet"])`

- 두 개 이상의 그림 비교해서 보기: `plot_comp(list_of_sd[, plot="image", figsize=default, perc=100, cmap="gray_r", fill=True, scale=1, dx=0, key=False])`

- 두 개의 데이터 트레이스 비교해서 그리기(tnum: trace number (-1: 전체 평균 트레이스): `seis_env_dB(sd, gained_sd[, tnum=-1])`

- 송수신기 좌표 보기: `stacking_chart(sd)`
- Autocorrelation map: `auto_correlation_map(sd[, double max_lag=0.2])`

- 속도모델 그리기: `plot_vel(vel, h[, figsize=[15,4], unit="km/s", tick=np.arange(1.5,6,1)])`
- 구조보정 결과 그리기: `plot_mig(mig, h[, figsize=[15,4]])`

## Processing

- key 값을 기준으로 정렬(list_of_keywords에서 부호: (+)오름차순, (-)내림차순  ex. ["+fldr","-offset"])(susort): `sorted_sd = trace_sort(sd, list_of_keywords)`

- 이득 조절 적용하기(sugain): `gained_sd = gain(sd[, tpow=0, epow=0, agc=False, agc_gate=0.5, norm="rms"|"amplitude"])`

- 주파수 필터링(cut_off=[min_freq, max_freq])(sufilter): `filtered_sd = bpfilter(sd, cut_off)`

- 자기상관: `corr_trc = autocorr(double trc[:])`
- 송곳 곱풀기: `ddata = spiking_decon(sd[, double max_lag=0.2, double mu=0.1])`
- 정적 보정(Surface-consistent static correction): `sdata = scr_static(sd, int cmp_start, int cmp_end, double maxlags)`

- 겹쌓기: `stacked_sd = stack(sd)`
- Stolt 참반사 보정(구조보정, v=상속도, dx=트레이스 간격): `mig_sd = stolt_mig(sd, v, dx)`

- 주시 계산: `time = traveltime(double vel[nx,nz], double h, int srcx, int srcz)`
- Kirchhoff prestack migration: `image = kirchhoff(sd, h, times[nx,nx,nz], tdelay)`
- 깊이 방향 2차 미분: `dimage = zdiff2(image)`
- 속도모델 평활화: `svel = moving_average2d(vel, r1, r2)`

## Velocity Analysis

### GUI
- Left panel: CMP gather before NMO
- Middle panel: Semblance panel
- Right panel: CMP gather after NMO

### Picking on the semblance panel
- 마우스 왼쪽 버튼 클릭: 점 찍기
- 마우스 오른쪽 버튼 클릭: 점 지우기

## 탄성파 자료 직접 접근
- 트레이스 헤더: sd.header
- 시계열 자료: sd.data
