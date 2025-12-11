'use client';

import { useEffect, useRef, useState } from 'react';
import { DisplayData } from '@/types/terminal';
import { displayMatrix, displayColorMatrix } from '@/lib/matrix-utils';

declare global {
  interface Window {
    MathJax: any;
  }
}

interface MatrixOutputProps {
  data?: DisplayData | null;
}

export default function MatrixOutput({ data }: MatrixOutputProps) {
  const outputRef = useRef<HTMLDivElement>(null);
  const [outputs, setOutputs] = useState<string[]>([]);

  useEffect(() => {
    // Load MathJax
    if (typeof window !== 'undefined' && !window.MathJax) {
      const script = document.createElement('script');
      script.src = 'https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML';
      script.async = true;
      document.head.appendChild(script);

      window.MathJax = {
        tex: {
          inlineMath: [['$', '$'], ['\\(', '\\)']],
          displayMath: [['$$', '$$'], ['\\[', '\\]']],
        },
        svg: {
          fontCache: 'global',
        },
      };
    }
  }, []);

  useEffect(() => {
    if (!data) return;

    let latex = '';
    const varName = 'var'; // You can pass this as a prop

    switch (data.what) {
      case 'matrix':
        if (data.value) {
          latex = `$$${displayMatrix(data.value)}$$`;
        }
        break;

      case 'color_matrix':
        if (data.value) {
          latex = `$$${displayColorMatrix(data.value)}$$`;
        }
        break;

      case 'ss':
        if (data.A && data.B && data.C && data.D) {
          latex = `$$ \\dot{x} = ${displayMatrix(data.A.value)}x + ${displayMatrix(data.B.value)}u $$`;
          latex += `$$ y = ${displayMatrix(data.C.value)}x + ${displayMatrix(data.D.value)}u $$`;
        }
        break;

      case 'orbit':
        latex = `\\begin{eqnarray} h & = & ${data.h}\\,km^2/s \\nonumber \\\\`;
        latex += `a & = & ${data.a}\\, km \\nonumber \\\\`;
        latex += `e & = & ${data.e}\\nonumber \\\\`;
        latex += `\\Omega & = & ${data.Omega}^{\\circ} \\nonumber \\\\`;
        latex += `i & = & ${data.i}^{\\circ} \\nonumber \\\\`;
        latex += `\\omega & = & ${data.omega}^{\\circ} \\nonumber \\end{eqnarray}`;
        break;

      case 'controller':
        if (data.A && data.B1 && data.B2 && data.C && data.D1 && data.D2) {
          latex = `$$ \\dot{x}_c = ${displayMatrix(data.A.value)}x + ${displayMatrix(data.B1.value)}y + ${displayMatrix(data.B2.value)}r$$`;
          latex += `$$ u = ${displayMatrix(data.C.value)}x_c + ${displayMatrix(data.D1.value)}y + ${displayMatrix(data.D2.value)}r$$`;
        }
        break;

      case 'feedback':
        if (data.A && data.B && data.C && data.D) {
          latex = `$$ \\dot{x}_a = ${displayMatrix(data.A.value)}x_a + ${displayMatrix(data.B.value)}r $$`;
          latex += `$$ y = ${displayMatrix(data.C.value)}x_a + ${displayMatrix(data.D.value)}r $$`;
        }
        break;

      case 'plot':
        if (data.script && data.div) {
          // Handle Bokeh plots
          const plotHtml = data.script + '<br />' + data.div;
          setOutputs(prev => [plotHtml, ...prev]);
          return;
        }
        break;

      case 'image':
        if (data.value) {
          const imgHtml = `<img src="/uploads/${data.value}" width="500" alt="Uploaded image" class="rounded-lg shadow-lg"/>`;
          setOutputs(prev => [imgHtml, ...prev]);
          return;
        }
        break;
    }

    if (latex) {
      setOutputs(prev => [latex, ...prev]);

      // Trigger MathJax rendering
      if (window.MathJax && window.MathJax.typesetPromise) {
        setTimeout(() => {
          window.MathJax.typesetPromise([outputRef.current]);
        }, 100);
      } else if (window.MathJax && window.MathJax.Hub) {
        setTimeout(() => {
          window.MathJax.Hub.Queue(['Typeset', window.MathJax.Hub, outputRef.current]);
        }, 100);
      }
    }
  }, [data]);

  return (
    <div className="mt-8">
      <h2 className="text-2xl font-bold mb-4 text-gray-200">Graphical Output</h2>
      <div
        ref={outputRef}
        className="bg-gray-800 border border-gray-700 rounded-lg p-6 min-h-[200px] space-y-4"
      >
        {outputs.length === 0 ? (
          <p className="text-gray-500 italic">
            No output yet. Use gdisplay command to show matrices and plots.
          </p>
        ) : (
          outputs.map((output, index) => (
            <div
              key={index}
              className="border-b border-gray-700 last:border-0 pb-4 last:pb-0"
              dangerouslySetInnerHTML={{ __html: output }}
            />
          ))
        )}
      </div>
    </div>
  );
}
