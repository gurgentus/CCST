# FastAPI router for control operations
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional
import control  # Import the pure Python control module

router = APIRouter()

# Request models
class MatrixRequest(BaseModel):
    name: str
    value: str  # JSON string of matrix

class ModesRequest(BaseModel):
    mat: str
    V: str
    E: str

class SSRequest(BaseModel):
    A: str
    B: str
    C: str
    D: str
    name: str

class ControllerRequest(BaseModel):
    A: str
    B1: str
    B2: str
    C: str
    D1: str
    D2: str
    name: str

class FeedbackRequest(BaseModel):
    G: str
    K: str
    name: str

class ControlRequest(BaseModel):
    index: int
    G: str
    name: str

class OutputRequest(BaseModel):
    index: int
    G: str
    name: str

class ResponseRequest(BaseModel):
    output: str
    G: str
    init: str
    inp: str
    time: float
    plot: str

class StepRequest(BaseModel):
    output: str
    G: str
    control: str
    st: float
    time: float
    plot: str


# Endpoints - thin layer that calls control module functions
@router.post("/matrix")
async def matrix_endpoint(req: MatrixRequest):
    """Define a matrix"""
    try:
        result = control.matrix(req.name, req.value)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/modes")
async def modes_endpoint(req: ModesRequest):
    """Calculate eigenvectors and eigenvalues of a matrix"""
    try:
        result = control.modes(req.mat, req.V, req.E)
        return {"result": result}
    except KeyError as e:
        raise HTTPException(status_code=404, detail=f"Matrix '{req.mat}' not found")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/ss")
async def ss_endpoint(req: SSRequest):
    """Define state-space system"""
    try:
        result = control.ss(req.A, req.B, req.C, req.D, req.name)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/controller")
async def controller_endpoint(req: ControllerRequest):
    """Define a controller"""
    try:
        result = control.controller(req.A, req.B1, req.B2, req.C, req.D1, req.D2, req.name)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/feedback")
async def feedback_endpoint(req: FeedbackRequest):
    """Create feedback system"""
    try:
        result = control.feedback(req.G, req.K, req.name)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/control")
async def control_endpoint(req: ControlRequest):
    """Set control variable"""
    try:
        result = control.control(req.index, req.G, req.name)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/output")
async def output_endpoint(req: OutputRequest):
    """Set output variable"""
    try:
        result = control.output(req.index, req.G, req.name)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/response")
async def response_endpoint(req: ResponseRequest):
    """Calculate system response"""
    try:
        result = control.response(req.output, req.G, req.init, req.inp, req.time, req.plot)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/step")
async def step_endpoint(req: StepRequest):
    """Calculate step response"""
    try:
        result = control.step(req.output, req.G, req.control, req.st, req.time, req.plot)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
