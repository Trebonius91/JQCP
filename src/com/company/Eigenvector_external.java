package com.company;

public class Eigenvector_external {
    /**
	 ************
	 * Array transpose.
	 *
	 * @param paraArray The transposing array(double[][]).
	 *
	 * @return returnArray The transposed array(double[][]).
	 ************
	 */
	public static double[][] arrayTranspose(double[][] paraArray) {
		int tempM = paraArray.length;
		int tempN = paraArray[0].length;
		double[][] returnArray = new double[tempN][tempM];
		for (int i = 0; i < tempM; i++) {
			for (int j = 0; j < tempN; j++) {
				returnArray[j][i] = paraArray[i][j];
			} // Of for j
		} // Of for i

		return returnArray;
	}// Of arrayTranspose
    /**
	 ************
	  * Array corresponding elements multiply and add all elements.
	 *
	 * @param paraFirstArray  The first array(double[]).
	 * @param paraSecondArray The second array(double[]).
	 * @return The result(double).
	 ************
	 */
	public static double arrayMultiplyAndAdd(double[] paraFirstArray, double[] paraSecondArray) {
		int tempM = paraFirstArray.length;
		double resultMultipliedArray = 0;

		for (int i = 0; i < tempM; i++) {
			resultMultipliedArray += paraFirstArray[i] * paraSecondArray[i];
		} // Of for i

		return resultMultipliedArray;
	}// Of arrayMultiplyAndAdd

	/**
	 ************
	 * Get the index of array.
	 *
	 * @param paraLen The length of array(int).
	 * @return The index of array(int[]).
	 ************
	 */
	public static int[] arrayIndexAuto(int paraLen) {
		int[] returnArray = new int[paraLen];
		for (int i = 0; i < paraLen; i++) {
			returnArray[i] = i;
		}
		return returnArray;
	}// Of arrayIndexAuto

	/**
	 ************
	 * Get the index of array.
	 *
	 * @param paraStrat The start index(int).
	 * @param paraEnd   The end index(int).
	 * @return The index of array(int[]).
	 ************
	 */
	public static int[] arrayIndexAuto(int paraStrat, int paraEnd) {
		int[] returnArray = new int[paraEnd - paraStrat];
		for (int i = 0; i < returnArray.length; i++) {
			returnArray[i] = i + paraStrat;
		}
		return returnArray;
	}// Of arrayIndexAuto

	/**
	 ************
	 * Copy array.
	 *
	 * @param paraArray The array(double[][]).
	 * @return The copied array(double[][]).
	 ************
	 */
	public static double[][] arrayCopy(double[][] paraArray) {
		int tempM = paraArray.length;
		double[][] resultArray = new double[tempM][];
		for (int i = 0; i < tempM; i++) {
			resultArray[i] = paraArray[i];
		} // Of for i
		return resultArray;
	}// Of arrayCopy
    /**
	 ************
	 * Get the elements of the corresponding row of the incoming array.
	 *
	 * @param paraArray The import array(double[][]).
	 * @param paraIndex The index i(int[]).
	 * @return The row elements corresponding incoming array(doubel[][]).
	 ************
	 */
	public static double[][] arrayRowValue(double[][] paraArray, int[] paraIndex) {
		int tempParaIndex = paraIndex.length;
		double[][] returnRowValue = new double[tempParaIndex][];

		for (int i = 0; i < returnRowValue.length; i++) {
			returnRowValue[i] = paraArray[paraIndex[i]];
		} // Of for i

		return returnRowValue;
	}// Of arrayRowValue

	/**
	 ************
	  * The determinant of matrix. Solve matrix determinant
	 *
	 * @param paraMatrix The given matrix(double[][]).
	 *
	 * @return The answer of determinant.
	 ************
	 */

	public static double matrixDeterminant(double[][] paraMatrix) {

		Matrix matrix = new Matrix();
		int tempM = paraMatrix.length;
		int tempN = paraMatrix[0].length;

		System.out.println("The para matrix");
		matrix.print(paraMatrix);
		if (tempM == 1) {
			return paraMatrix[0][0];
		} else if (tempM == 2) {
			return paraMatrix[0][0] * paraMatrix[1][1] - paraMatrix[1][0] * paraMatrix[0][1];
		} // Of if

		double resultValue = 0;
		int k;
		for (int i = 0; i < tempN; i++) {
			double[][] tempMatrix = new double[tempM - 1][tempN - 1];
			for (int j = 0; j < tempM - 1; j++) {
				k = 0;
				while (k < tempN - 1) {
					if (k < i) {
						tempMatrix[j][k] = paraMatrix[j + 1][k];
					} else {
						tempMatrix[j][k] = paraMatrix[j + 1][k + 1];
					} // Of if
					k++;
				} // Of while
			} // Of for j
			resultValue += (Math.pow(-1., i)) * paraMatrix[0][i] * matrixDeterminant(tempMatrix);
		} // Of for i

		return resultValue;
	}// Of matrixDeterminant

    /**
	 ************
	  * The gram schimidt of matrix. Solve matrix orthogonalization
	 *
	 * @param paraMatrix The given matrix(double[][]).
	 *
	 * @return The gram schimidt of matrix.
	 ************
	 */
	public static double[][] matrixGramSchimidt(double[][] paraMatrix) {
		if (paraMatrix == null) {
			return null;
		} // Of if

		double[][] tempTransposedMatrix = arrayTranspose(paraMatrix);
		int tempM = tempTransposedMatrix.length;
		int tempN = tempTransposedMatrix[0].length;

		double[][] resultMatrix = new double[tempM][tempN];
		double tempValue = 0;
		double tempFactor = 0;
		for (int i = 0; i < tempM; i++) {
			for (int j = 0; j < tempN; j++) {
				tempValue = tempTransposedMatrix[i][j];
				for (int k = 0; k < i; k++) {
					tempFactor = (1. * arrayMultiplyAndAdd(tempTransposedMatrix[i], resultMatrix[k]))
							/ arrayMultiplyAndAdd(resultMatrix[k], resultMatrix[k]);
					tempValue -= tempFactor * resultMatrix[k][j];
				} // Of for k
				resultMatrix[i][j] = tempValue;
			} // Of for j
		} // Of for i

		return arrayTranspose(resultMatrix);
	}// Of matrixGramSchimidt

    /**
	 ************
	  * The inverse of matrix.
	 *
	 * @param paraMatrix The given matrix(double[][]).
	 *
	 * @return The inverse of matrix(double[][]).
	 ************
	 */
	public static double[][] matrixInverse(double[][] paraMatrix) {
		if (paraMatrix == null) {
			return null;
		} // Of if

		// Make sure that this matrix is invertible.
		if (matrixDeterminant(paraMatrix) == 0) {
			System.out.println("Fetal error, the matrix can not be invertible.");
			System.exit(0);
		} // Of if

		int tempM = paraMatrix.length;
		int tempN = paraMatrix[0].length * 2;

		double[][] tempCopyMatrix = new double[tempM][tempN];

		for (int i = 0; i < tempM; i++) {
			for (int j = 0; j < tempN; j++) {
				if (j < tempM) {
					tempCopyMatrix[i][j] = paraMatrix[i][j];
				} else {
					if (j == i + tempM) {
						tempCopyMatrix[i][j] = 1;
					} // Of if
				} // Of if
			} // Of for j
		} // Of for i
       double tempTimes = 0;
		for (int i = 0; i < tempM; i++) {
			if (tempCopyMatrix[i][i] == 0) {
				int j;
				for (j = i + 1; j < tempM; j++) {
					if (tempCopyMatrix[j][i] != 0) {
						break;
					} // Of if
				} // Of for j

				if (j != i + 1) {
					for (int k = 0; k < tempN; k++) {
						double temppVlaue = tempCopyMatrix[i][k];
						tempCopyMatrix[i][k] = tempCopyMatrix[j][k];
						tempCopyMatrix[j][k] = temppVlaue;
					} // Of for k
				} // Of if
			} // Of if

			for (int j = i + 1; j < tempM; j++) {
				if (tempCopyMatrix[j][i] != 0) {
					tempTimes = (tempCopyMatrix[j][i]) / tempCopyMatrix[i][i];
					for (int k = i; k < tempN; k++) {
						tempCopyMatrix[j][k] /= tempTimes;
						tempCopyMatrix[j][k] -= tempCopyMatrix[i][k];
					} // Of for k
				} // Of if
			} // Of for j
		} // Of for i
        for (int i = 0; i < tempM; i++) {
			for (int j = i + 1; j < tempN / 2.; j++) {
				if (tempCopyMatrix[i][j] != 0) {
					tempTimes = tempCopyMatrix[i][j] / tempCopyMatrix[j][j];
					for (int k = j; k < tempN; k++) {
						tempCopyMatrix[i][k] -= tempTimes * tempCopyMatrix[j][k];
					} // Of for k
				} // Of if
			} // Of for j
		} // Of for i

		for (int i = 0; i < tempM; i++) {
			tempTimes = tempCopyMatrix[i][i];
			for (int j = 0; j < tempN; j++) {
				tempCopyMatrix[i][j] /= tempTimes;
			} // Of for j
		} // Of for i

		double[][] resultMatrix = new double[tempM][tempN / 2];
		for (int i = 0; i < tempM; i++) {
			for (int j = 0; j < tempN / 2; j++) {
				resultMatrix[i][j] = tempCopyMatrix[i][j + tempN / 2];
			} // Of for j
		} // Of for i

		return resultMatrix;
	}// Of matrixInverse
    /**
	 ************
	  * The QR-decomposition of matrix. Matrix QR decomposition
	 *
	 * @param paraMatrix  The given matrix(double[][]).
	 //* @param paraMatrixQ The given matrix, the size echo paraMatrix(double[][]).
	 ************
	 */
	public static double[][] matrixQrDecomposition(double[][] paraMatrix) {
		double[][] tempOrthogonalMatrix = arrayTranspose(matrixGramSchimidt(paraMatrix));
		int tempM = tempOrthogonalMatrix.length;
		int tempN = tempOrthogonalMatrix[0].length;

		double[][] tempMatrixQ = new double[tempM][tempN];
		for (int i = 0; i < tempM; i++) {
			double tempMag = magnitude(tempOrthogonalMatrix[i]);
			for (int j = 0; j < tempN; j++) {
				tempMatrixQ[i][j] = tempOrthogonalMatrix[i][j] / tempMag;
			} // Of for j
		} // Of for i

		double[][] tempMatrixR = matrixMultiply(tempMatrixQ, paraMatrix);
		double[][] resultSummary = new double[tempM + tempN][tempM];
		for (int i = 0; i < tempN; i++) {
			for (int j = 0; j < tempM; j++) {
				resultSummary[i][j] = tempMatrixQ[j][i];
			} // Of for j
		} // Of for i

		for (int i = tempN; i < resultSummary.length; i++) {
			for (int j = 0; j < tempM; j++) {
				resultSummary[i][j] = tempMatrixR[i - tempN][j];
			} // Of for j
		} // Of for i

		return resultSummary;
	}// Of matrixQrDecomposition
    /**
	 ************
	 * The magnitude.
	 *
	 * @param paraMatrix The given matrix(double[]).
	 * @return The answer.
	 ************
	 */
	private static double magnitude(double[] paraMatrix) {
		return Math.sqrt(arrayMultiplyAndAdd(paraMatrix, paraMatrix));
	}// Of magnitude

	/**
	 ************
	 * The matrix multiply.
	 *
	 * @param paraFirstMatrix  The given first matrix(double[][]).
	 * @param paraSecondMatrix The given second matrix(double[][]).
	 * @return The answer(double[][]).
	 ************
	 */
	public static double[][] matrixMultiply(double[][] paraFirstMatrix, double[][] paraSecondMatrix) {
		if (paraFirstMatrix == null || paraSecondMatrix == null) {
			return null;
		} // Of if

		if (paraFirstMatrix[0].length != paraSecondMatrix.length) {
			System.out.println("Fetal error: The size of two inputed matrix is illegally.");
			System.exit(0);
		} // Of if

		int tempMa = paraFirstMatrix.length;
		int tempNa = paraFirstMatrix[0].length;
		int tempNb = paraSecondMatrix[0].length;
		double[][] resultMatrix = new double[tempMa][tempNb];
		for (int i = 0; i < tempMa; i++) {
			for (int j = 0; j < tempNb; j++) {
				double tempSum = 0;
				for (int k = 0; k < tempNa; k++) {
					tempSum += paraFirstMatrix[i][k] * paraSecondMatrix[k][j];
				} // Of for k
				resultMatrix[i][j] = tempSum;
			} // Of for j
		} // Of for i
		return resultMatrix;
	}// Of matrixMultiply

	/**
	 ************
	 * The eig of matrix.
	 *
	 * @param paraMatrix The given matrix(double[][]).
	 * @param paraIter   The maximum iteration.
	 * @return The eig of matrix.
	 ************
	 */
	public static double[][] matrixEigValue(double[][] paraMatrix, int paraIter) {
		int tempM = paraMatrix.length;
		int tempN = paraMatrix[0].length;
		int[] tempIndexQ = arrayIndexAuto(tempM);
		int[] tempIndexR = arrayIndexAuto(tempM, tempM + tempN);
		for (int i = 0; i < paraIter; i++) {
			double[][] tempSummary = matrixQrDecomposition(paraMatrix);
			double[][] tempMatrixQ = arrayRowValue(tempSummary, tempIndexQ);
			double[][] tempMatrixR = arrayRowValue(tempSummary, tempIndexR);
			paraMatrix = matrixMultiply(tempMatrixR, tempMatrixQ);
		} // Of for i

		return paraMatrix;
	}// Of matrixEigValue

	/**
	 ************
	  * The eig vector of matrix.
	 *
	 * @param paraMatrix The given matrix(double[][]).
	 * @param paraIter   The maximum iteration.
	 * @return The eig vector of matrix.
	 ************
	 */
	public static double[][] matrixEigVector(double[][] paraMatrix, int paraIter) {
		int tempM = paraMatrix.length;
		int tempN = paraMatrix[0].length;
		double[][] tempMatrix = matrixEigValue(arrayCopy(paraMatrix), paraIter);
		for (int i = 0; i < tempM; i++) {
			for (int j = 0; j < tempN; j++) {
				if (i != j) {
					tempMatrix[i][j] = 0;
				}//Of if
			}//Of for j
		}//Of for i

		return matrixInverse(matrixSubtract(paraMatrix, tempMatrix));
	}//Of matrixEigVector

	/**
	 ************
	  * The subtract of matrix.
	 *
	 * @param paraFirstMatrix The first given matrix(double[][]).
	 * @param paraFirstMatrix The second given matrix(double[][]).
	 * @return The subtract of two matrix.
	 ************
	 */
	public static double[][] matrixSubtract(double[][] paraFirstMatrix, double[][] paraSecondMatrix) {
		if (paraFirstMatrix == null || paraSecondMatrix == null) {
			return null;
		}//Of if

		int tempMFirst = paraFirstMatrix.length;
		int tempNFirst = paraFirstMatrix[0].length;
		int tempMSecond = paraSecondMatrix.length;
		int tempNSecond = paraSecondMatrix[0].length;
		if (tempMFirst != tempMSecond || tempNFirst != tempNSecond) {
			System.out.println("Fetal error, the inputed matrix is illegal in matrixSub(double[][], double[][])");
			System.exit(0);
		}//Of if

		double[][] resultMatrix = new double[tempMFirst][tempNFirst];
		for (int i = 0; i < tempMFirst; i++) {
			for (int j = 0; j < tempNFirst; j++) {
				resultMatrix[i][j] = paraFirstMatrix[i][j] - paraSecondMatrix[i][j];
			}//Of for j
		}//Of for i

		return resultMatrix;
	}//Of matrixSubtract

}
