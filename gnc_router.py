# FastAPI router for GNC operations
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
import helper

router = APIRouter()

# TODO: Port gnc.py methods to FastAPI
# For now, just placeholder endpoints

@router.get("/")
async def gnc_index():
    """GNC API index"""
    return {"message": "GNC API - methods to be implemented"}
