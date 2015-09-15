# IPythonProcessing

본 패키지는 2차원 탄성파 자료 처리 학부 수업을 위한 패키지입니다.

## 설치

1. 먼저 Python(version 2.7)과 관련 패키지들을 설치하거나(참고: [Anaconda](http://continuum.io/downloads)) Python을 사용할 수 있는 웹사이트 계정을 만드세요(참고: [SageMathCloud](https://cloud.sagemath.com), [Wakari](https://www.wakari.io)).
	1. 필요한 Python 모듈: Python, IPython, Numpy, SciPy, Matplotlib, Cython (Anaconda를 이용하면 한 번에 설치 가능)
	2. 주의: 한글 윈도우에서 사용자 이름이 한글일 때 IPython notebook 실행에 문제가 발생한 적 있음.
	3. 주의: SageMathCloud나 Wakari 무료 계정에서는 파일 크기 제한으로 인해 Marmousi 자료 전체(69MB)를 한 번에 처리할 수 없음.
2. [Github](https://github.com/pkgpl/IPythonProcessing)에서 "Download ZIP" 링크를 통해 패키지를 받고 원하는 위치(디렉토리 또는 웹사이트)에 압축을 푸세요.
3. pkprocess 디렉토리에서 `$python setup_cy.py build_ext --inplace` 명령을 실행해서 Cython 모듈을 설치하세요.
4. [IPython Notebook](http://ipython.org/notebook.html)을 이용해 압축을 푼 디렉토리에서 예제 파일을 열어 직접 실행해보세요.
	1. pkprocess 디렉토리를 PYTHON_PATH 환경변수에 추가하면 임의의 위치에서 패키지를 import할 수 있습니다.


## 예제

1. [Land data](http://nbviewer.ipython.org/github/pkgpl/IPythonProcessing/blob/master/Example_Land.ipynb) (Notebook 용량: 18 MB)
2. [Marine data](http://nbviewer.ipython.org/github/pkgpl/IPythonProcessing/blob/master/Example_Marine.ipynb) (Marmousi data, Notebook 용량: 10 MB)
3. [Marine data](http://nbviewer.ipython.org/github/pkgpl/IPythonProcessing/blob/master/Kirchhoff_mig.ipynb) (Kirchhoff migration, Notebook 용량: 1 MB 이하)

* 예제 문서는 IPython Notebook version 3.0 이상을 기준으로 만들었습니다. 예제에 사용한 데이터 파일(Land data: 3.5 MB, Marmousi data: 69 MB, Marmousi velocity model: 1 MB 이하)은 여기에서 제공하지 않습니다. 필요하신 분은 메일 주세요.


## Python 사용법

1. IPython Notebook [홈페이지](http://ipython.org/notebook.html), [동영상](https://www.youtube.com/watch?v=lZChNqEYqLA)
2. Python [홈페이지](https://www.python.org), [소개](https://www.wakari.io/nb/url///wakari.io/static/notebooks/Lecture_1_Introduction_to_Python_Programming.ipynb)
3. Numpy [홈페이지](http://www.numpy.org), [소개](https://www.wakari.io/nb/url///wakari.io/static/notebooks/Lecture_2_Numpy.ipynb)
4. Matplotlib [홈페이지](http://matplotlib.org), [소개](https://www.wakari.io/nb/url///wakari.io/static/notebooks/Lecture_4_Matplotlib.ipynb)


## 자료 처리 패키지 명령어

[매뉴얼](./Manual.md)을 참고하세요.


## Reference
하완수, 2015, 대화식 탄성파 자료 처리 수업을 위한 파이썬 패키지 개발, 한국 자원공학회지, Vol. 52, No. 4, pp. 414-421.
