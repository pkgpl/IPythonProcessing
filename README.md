# IPythonProcessing

본 패키지는 한국자원공학회지 "대화식 탄성파 자료 처리 수업을 위한 파이썬 패키지 개발" 논문에서 소개한 학부 2차원 탄성파 자료 처리 수업을 위한 패키지입니다.

## 설치
1. 먼저 Python(version 2.7)과 관련 패키지들을 설치하거나(참고: [Anaconda](http://continuum.io/downloads)) Python을 사용할 수 있는 웹사이트 계정을 만드세요(참고: [SageMathCloud](https://cloud.sagemath.com), [Wakari](https://www.wakari.io)).
2. [Github](https://github.com/pkgpl/IPythonProcessing)에서 "Download ZIP" 링크를 통해 패키지를 받고 원하는 위치(디렉토리 또는 웹사이트)에 압축을 푸세요.
3. pkprocess 디렉토리에서 `$python setup_cy.py build_ext --inplace` 명령을 실행해서 Cython 모듈을 설치하세요.

## 예제
1. [Marine data](https://github.com/pkgpl/IPythonProcessing/blob/master/Example_Marine.ipynb) (Marmousi data)
2. [Land data](https://github.com/pkgpl/IPythonProcessing/blob/master/Example_Land.ipynb)

* 예제에 사용한 데이터 파일은 여기에서 제공하지 않습니다. 필요하신 분은 메일 주세요.