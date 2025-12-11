// Matrix display utilities for LaTeX rendering

export function displayMatrix(matrix: number[][]): string {
  let latex = "\\begin{pmatrix}";

  for (let i = 0; i < matrix.length; i++) {
    for (let j = 0; j < matrix[i].length; j++) {
      latex += matrix[i][j];
      if (j < matrix[i].length - 1) {
        latex += "&";
      }
    }
    if (i < matrix.length - 1) {
      latex += "\\\\";
    }
  }

  latex += "\\end{pmatrix}";
  return latex;
}

export function displayColorMatrix(matrix: number[][]): string {
  let latex = "\\begin{pmatrix}";

  for (let i = 0; i < matrix.length; i++) {
    for (let j = 0; j < matrix[i].length; j++) {
      let color = "{";
      if (matrix[i][j] === 2) {
        color = "{\\colorbox{green}";
      }
      if (matrix[i][j] === 1) {
        color = "{\\colorbox{red}";
      }
      latex += color + matrix[i][j] + "}";
      if (j < matrix[i].length - 1) {
        latex += "&";
      }
    }
    if (i < matrix.length - 1) {
      latex += "\\\\";
    }
  }

  latex += "\\end{pmatrix}";
  return latex;
}

export function echoMatrix(matrix: number[][]): string {
  let text = "";

  for (let i = 0; i < matrix.length; i++) {
    for (let j = 0; j < matrix[i].length; j++) {
      text += matrix[i][j];
      if (j < matrix[i].length - 1) {
        text += " ";
      }
    }
    if (i < matrix.length - 1) {
      text += "\n";
    }
  }

  return text;
}
